#include "mds/mining/orbit.hpp"

#include <algorithm>
#include <bit>
#include <limits>
#include <vector>

#include "mds/flat_bitmap.hpp"
#include "mds/fxhash.hpp"
#include "mds/mining/core.hpp"
#include "mds/profiler.hpp"
#include "mds/simd_bitmap.hpp"

namespace mds {

namespace {

template <typename MarkPredicate>
struct OrbitBitmapInfo {
    FlatBitmap bitmaps;
    std::vector<std::size_t> orbit_size;
    std::vector<std::int32_t> orbit_label;
    std::vector<std::size_t> orbit_count;
    std::vector<std::vector<std::size_t>> label_to_orbits;
    std::size_t bitmap_words = 0;
};

std::size_t find_root(std::vector<std::size_t>& parent, std::size_t x) {
    if (parent[x] != x) {
        parent[x] = find_root(parent, parent[x]);
    }
    return parent[x];
}

void union_root(std::vector<std::size_t>& parent, std::size_t x, std::size_t y) {
    const auto px = find_root(parent, x);
    const auto py = find_root(parent, y);
    if (px != py) {
        parent[px] = py;
    }
}

std::vector<std::vector<std::size_t>> find_automorphisms(const Graph& q) {
    const auto n = static_cast<std::size_t>(q.n_vertex);
    std::vector<std::vector<std::size_t>> result;

    auto max_label = 0;
    for (std::size_t u = 0; u < n; ++u) {
        max_label = std::max(max_label, q.label[u]);
    }
    std::vector<std::vector<std::size_t>> label_groups(static_cast<std::size_t>(max_label + 1));
    for (std::size_t u = 0; u < n; ++u) {
        label_groups[static_cast<std::size_t>(q.label[u])].push_back(u);
    }

    std::vector<std::int32_t> adj_matrix(n * n, -1);
    for (std::size_t u = 0; u < n; ++u) {
        for (auto i = static_cast<std::size_t>(q.nbr_offset[u]); i < static_cast<std::size_t>(q.nbr_offset[u + 1]); ++i) {
            const auto v = static_cast<std::size_t>(q.nbr[i].first);
            adj_matrix[u * n + v] = q.nbr[i].second;
        }
    }

    std::vector<std::size_t> perm(n, static_cast<std::size_t>(-1));
    std::vector<bool> used(n, false);

    const auto backtrack = [&](auto&& self, std::size_t u) -> void {
        if (u == n) {
            result.push_back(perm);
            return;
        }

        const auto label = q.label[u];
        const auto label_idx = static_cast<std::size_t>(label);
        if (label_idx >= label_groups.size()) {
            return;
        }

        for (const auto v : label_groups[label_idx]) {
            if (used[v]) {
                continue;
            }

            bool valid = true;
            for (std::size_t prev_u = 0; prev_u < u; ++prev_u) {
                const auto prev_v = perm[prev_u];
                if (adj_matrix[prev_u * n + u] != adj_matrix[prev_v * n + v] ||
                    adj_matrix[u * n + prev_u] != adj_matrix[v * n + prev_v]) {
                    valid = false;
                    break;
                }
            }

            if (valid) {
                perm[u] = v;
                used[v] = true;
                self(self, u + 1);
                used[v] = false;
            }
        }
    };

    backtrack(backtrack, 0);
    if (result.empty()) {
        result.emplace_back();
        result.back().reserve(n);
        for (std::size_t i = 0; i < n; ++i) {
            result.back().push_back(i);
        }
    }
    return result;
}

bool has_any_intersection(const std::uint64_t* lhs, const std::uint64_t* rhs, std::size_t words) {
    for (std::size_t i = 0; i < words; ++i) {
        if ((lhs[i] & rhs[i]) != 0ULL) {
            return true;
        }
    }
    return false;
}

std::vector<std::vector<std::size_t>> build_connected_components(const std::vector<std::uint8_t>& intersect,
                                                                 std::size_t m) {
    if (m <= 1) {
        return m == 1 ? std::vector<std::vector<std::size_t>>{{0}} : std::vector<std::vector<std::size_t>>{};
    }

    std::vector<std::uint8_t> visited(m, 0);
    std::vector<std::vector<std::size_t>> components;
    for (std::size_t start = 0; start < m; ++start) {
        if (visited[start]) {
            continue;
        }

        std::vector<std::size_t> component;
        std::vector<std::size_t> stack{start};
        visited[start] = 1;

        while (!stack.empty()) {
            const auto u = stack.back();
            stack.pop_back();
            component.push_back(u);
            for (std::size_t v = 0; v < m; ++v) {
                if (!visited[v] && intersect[u * m + v] != 0) {
                    visited[v] = 1;
                    stack.push_back(v);
                }
            }
        }

        if (!component.empty()) {
            components.push_back(std::move(component));
        }
    }

    return components;
}

template <typename MarkPredicate>
OrbitBitmapInfo<MarkPredicate> build_orbit_bitmap_info(const Graph& q,
                                                       const std::vector<CandidateSpace>& cs,
                                                       const std::vector<std::size_t>& orbits,
                                                       std::size_t data_size,
                                                       MarkPredicate&& include_candidate) {
    const auto n = static_cast<std::size_t>(q.n_vertex);
    OrbitBitmapInfo<MarkPredicate> info;
    if (n == 0 || orbits.empty()) {
        return info;
    }

    std::size_t num_orbits = 0;
    for (const auto orbit_id : orbits) {
        num_orbits = std::max(num_orbits, orbit_id + 1);
    }
    if (num_orbits == 0) {
        return info;
    }

    info.bitmap_words = (data_size + 63) / 64;
    info.bitmaps = FlatBitmap(num_orbits, info.bitmap_words);
    info.orbit_size.assign(num_orbits, 0);
    info.orbit_label.assign(num_orbits, 0);
    info.orbit_count.assign(num_orbits, 0);

    std::vector<std::size_t> representative(num_orbits, n);
    for (std::size_t u = 0; u < n; ++u) {
        const auto oid = orbits[u];
        ++info.orbit_size[oid];
        if (representative[oid] == n) {
            representative[oid] = u;
        }
    }

    std::int32_t max_label = 0;
    for (std::size_t oid = 0; oid < num_orbits; ++oid) {
        const auto u = representative[oid];
        if (u >= n || u >= cs.size()) {
            continue;
        }

        const auto label = q.label[u];
        info.orbit_label[oid] = label;
        max_label = std::max(max_label, label);

        for (std::size_t i = 0; i < static_cast<std::size_t>(cs[u].size); ++i) {
            if (!include_candidate(cs[u].marked[i])) {
                continue;
            }
            info.bitmaps.set_bit(oid, static_cast<std::size_t>(cs[u].candidates[i]));
            ++info.orbit_count[oid];
        }
    }

    info.label_to_orbits.assign(static_cast<std::size_t>(max_label + 1), {});
    for (std::size_t oid = 0; oid < num_orbits; ++oid) {
        const auto label = info.orbit_label[oid];
        if (label >= 0) {
            info.label_to_orbits[static_cast<std::size_t>(label)].push_back(oid);
        }
    }

    return info;
}

std::int32_t compute_component_mds_bitmap(const std::vector<std::size_t>& component,
                                          const std::vector<std::size_t>& orbit_ids,
                                          const std::vector<std::size_t>& orbit_sizes,
                                          const FlatBitmap& bitmaps,
                                          std::size_t bitmap_words,
                                          std::vector<std::uint64_t>& temp_union,
                                          std::int32_t stop_at) {
    const auto m = component.size();
    if (m == 0) {
        return std::numeric_limits<std::int32_t>::max();
    }
    if (m == 1) {
        const auto oid = orbit_ids[component[0]];
        const auto size = orbit_sizes[oid];
        if (size == 0) {
            return std::numeric_limits<std::int32_t>::max();
        }
        return static_cast<std::int32_t>(bitmaps.popcount(oid) / size);
    }

    const auto subset_count = static_cast<std::size_t>(1) << m;
    std::vector<std::size_t> denoms(subset_count, 0);
    for (std::size_t mask = 1; mask < subset_count; ++mask) {
        const auto bit = static_cast<std::size_t>(std::countr_zero(mask));
        denoms[mask] = denoms[mask & (mask - 1)] + orbit_sizes[orbit_ids[component[bit]]];
    }

    std::int32_t best_mds = std::numeric_limits<std::int32_t>::max();
    for (std::size_t mask = 1; mask < subset_count; ++mask) {
        const auto total_size = denoms[mask];
        if (total_size == 0) {
            continue;
        }

        simd_bitmap::clear(temp_union.data(), bitmap_words);
        auto subset = mask;
        while (subset != 0) {
            const auto bit = static_cast<std::size_t>(std::countr_zero(subset));
            simd_bitmap::or_inplace(
                temp_union.data(),
                bitmaps.row(orbit_ids[component[bit]]),
                bitmap_words);
            subset &= (subset - 1);
        }

        const auto mds = static_cast<std::int32_t>(simd_bitmap::popcount(temp_union.data(), bitmap_words) / total_size);
        best_mds = std::min(best_mds, mds);
        if (best_mds <= stop_at) {
            return best_mds;
        }
    }

    return best_mds == std::numeric_limits<std::int32_t>::max() ? 0 : best_mds;
}

template <typename MarkPredicate>
std::int32_t compute_mds_with_bitmaps(const Graph& q,
                                      const std::vector<CandidateSpace>& cs,
                                      const std::vector<std::size_t>& orbits,
                                      std::size_t data_size,
                                      MarkPredicate&& include_candidate,
                                      std::optional<std::int32_t> stop_at) {
    const auto info = build_orbit_bitmap_info(q, cs, orbits, data_size, std::forward<MarkPredicate>(include_candidate));
    if (info.orbit_size.empty()) {
        return 0;
    }

    std::vector<std::uint64_t> temp_union(info.bitmap_words, 0ULL);
    std::int32_t min_mds = std::numeric_limits<std::int32_t>::max();
    const auto early_stop = stop_at.value_or(std::numeric_limits<std::int32_t>::min());

    for (const auto& orbit_ids : info.label_to_orbits) {
        const auto m = orbit_ids.size();
        if (m == 0) {
            continue;
        }

        if (m == 1) {
            const auto oid = orbit_ids[0];
            const auto size = info.orbit_size[oid];
            if (size > 0) {
                min_mds = std::min(min_mds, static_cast<std::int32_t>(info.orbit_count[oid] / size));
                if (stop_at.has_value() && min_mds <= early_stop) {
                    return min_mds;
                }
            }
            continue;
        }

        std::vector<std::uint8_t> intersect(m * m, 0);
        for (std::size_t i = 0; i < m; ++i) {
            intersect[i * m + i] = 1;
            for (std::size_t j = i + 1; j < m; ++j) {
                const auto has_intersection =
                    has_any_intersection(info.bitmaps.row(orbit_ids[i]), info.bitmaps.row(orbit_ids[j]), info.bitmap_words);
                const auto value = static_cast<std::uint8_t>(has_intersection ? 1 : 0);
                intersect[i * m + j] = value;
                intersect[j * m + i] = value;
            }
        }

        const auto components = build_connected_components(intersect, m);
        for (const auto& component : components) {
            const auto comp_mds = compute_component_mds_bitmap(component,
                                                               orbit_ids,
                                                               info.orbit_size,
                                                               info.bitmaps,
                                                               info.bitmap_words,
                                                               temp_union,
                                                               stop_at.value_or(std::numeric_limits<std::int32_t>::min()));
            min_mds = std::min(min_mds, comp_mds);
            if (stop_at.has_value() && min_mds <= early_stop) {
                return min_mds;
            }
        }
    }

    return min_mds == std::numeric_limits<std::int32_t>::max() ? 0 : min_mds;
}

}

std::vector<std::size_t> compute_orbits(const Graph& q) {
    const auto n = static_cast<std::size_t>(q.n_vertex);
    if (n == 0) {
        return {};
    }
    if (n == 1) {
        return {0};
    }

    std::vector<std::size_t> parent(n);
    for (std::size_t i = 0; i < n; ++i) {
        parent[i] = i;
    }

    const auto automorphisms = find_automorphisms(q);
    for (const auto& perm : automorphisms) {
        for (std::size_t u = 0; u < n; ++u) {
            union_root(parent, u, perm[u]);
        }
    }

    std::vector<std::size_t> root_to_orbit(n, n);
    std::vector<std::size_t> orbits(n, 0);
    std::size_t next_orbit = 0;
    for (std::size_t u = 0; u < n; ++u) {
        const auto root = find_root(parent, u);
        if (root_to_orbit[root] == n) {
            root_to_orbit[root] = next_orbit++;
        }
        orbits[u] = root_to_orbit[root];
    }

    return orbits;
}

std::int32_t compute_mds_ub(const Graph& q,
                            const std::vector<CandidateSpace>& cs,
                            const std::vector<std::size_t>& orbits,
                            MiningContext& ctx) {
    profile::ScopedPhase phase("mds_ub");
    return compute_mds_with_bitmaps(
        q,
        cs,
        orbits,
        static_cast<std::size_t>(ctx.data_graph.n_vertex),
        [](std::int8_t mark) { return mark >= 0; },
        ctx.tau);
}

std::int32_t compute_mds_from_cs(const Graph& q,
                                 const std::vector<CandidateSpace>& cs,
                                 const std::vector<std::size_t>& orbits,
                                 MiningContext& ctx) {
    return compute_mds_with_bitmaps(
        q,
        cs,
        orbits,
        static_cast<std::size_t>(ctx.data_graph.n_vertex),
        [](std::int8_t mark) { return mark == 1; },
        std::nullopt);
}

}
