#include "mds/mining/mds_context.hpp"

#include <algorithm>
#include <limits>

#include "mds/simd_bitmap.hpp"

namespace mds {

namespace bitmap {

inline void clear(std::vector<std::uint64_t>& bmp) {
    std::fill(bmp.begin(), bmp.end(), 0ULL);
}

inline void set_bit(std::vector<std::uint64_t>& bmp, std::size_t idx) {
    const auto word = idx / 64;
    const auto bit = idx % 64;
    if (word < bmp.size()) {
        bmp[word] |= (1ULL << bit);
    }
}

}

MDSContext::MDSContext() = default;

std::span<const std::size_t> MDSContext::orbit_members(std::size_t oid) const {
    if (oid >= num_orbits) {
        return {};
    }
    const auto start = orbit_members_offset_[oid];
    const auto size = orbit_size[oid];
    return {orbit_members_data_.data() + static_cast<std::ptrdiff_t>(start), size};
}

void MDSContext::init(const Graph& q,
                      std::vector<CandidateSpace>& cs,
                      const std::vector<std::size_t>& orbits_in,
                      std::size_t data_size) {
    bitmap_words_ = (data_size + 63) / 64;
    orbits = orbits_in;

    const auto max_orbit_it = std::max_element(orbits.begin(), orbits.end());
    num_orbits = max_orbit_it == orbits.end() ? 0 : (*max_orbit_it + 1);

    std::vector<std::size_t> orbit_sizes(num_orbits, 0);
    for (const auto orbit_id : orbits) {
        ++orbit_sizes[orbit_id];
    }

    orbit_members_offset_.assign(num_orbits, 0);
    std::size_t total_size = 0;
    for (std::size_t oid = 0; oid < num_orbits; ++oid) {
        orbit_members_offset_[oid] = total_size;
        total_size += orbit_sizes[oid];
    }

    orbit_members_data_.assign(total_size, 0);
    orbit_size = orbit_sizes;
    auto current_pos = orbit_members_offset_;
    for (std::size_t u = 0; u < orbits.size(); ++u) {
        orbit_members_data_[current_pos[orbits[u]]++] = u;
    }

    orbit_label.assign(num_orbits, 0);
    std::int32_t max_label = 0;
    for (std::size_t oid = 0; oid < num_orbits; ++oid) {
        const auto members = orbit_members(oid);
        if (!members.empty()) {
            const auto label = q.label[members[0]];
            orbit_label[oid] = label;
            max_label = std::max(max_label, label);
        }
    }

    num_labels = static_cast<std::size_t>(max_label + 1);
    label_to_orbits.assign(num_labels, {});
    for (std::size_t oid = 0; oid < num_orbits; ++oid) {
        const auto label = orbit_label[oid];
        if (label >= 0) {
            label_to_orbits[static_cast<std::size_t>(label)].push_back(oid);
        }
    }

    active_labels.clear();
    for (std::size_t label_idx = 0; label_idx < num_labels; ++label_idx) {
        if (!label_to_orbits[label_idx].empty()) {
            active_labels.push_back(label_idx);
        }
    }

    orbit_confirmed_ = FlatBitmap(num_orbits, bitmap_words_);
    orbit_remaining_ = FlatBitmap(num_orbits, bitmap_words_);
    orbit_confirmed_count_.assign(num_orbits, 0);
    orbit_remaining_count_.assign(num_orbits, 0);

    std::vector<std::uint64_t> temp_bitmap(bitmap_words_, 0ULL);
    for (std::size_t oid = 0; oid < num_orbits; ++oid) {
        const auto start = orbit_members_offset_[oid];
        const auto size = orbit_size[oid];
        if (size == 0) {
            continue;
        }

        bool first = true;
        for (std::size_t idx = 0; idx < size; ++idx) {
            const auto u = orbit_members_data_[start + idx];
            if (u >= cs.size()) {
                continue;
            }

            if (first) {
                for (std::size_t i = 0; i < static_cast<std::size_t>(cs[u].size); ++i) {
                    const auto v = static_cast<std::size_t>(cs[u].candidates[i]);
                    if (cs[u].marked[i] >= 0) {
                        orbit_remaining_.set_bit(oid, v);
                    }
                    if (cs[u].marked[i] == 1) {
                        orbit_confirmed_.set_bit(oid, v);
                    }
                }
                first = false;
            } else {
                bitmap::clear(temp_bitmap);
                for (std::size_t i = 0; i < static_cast<std::size_t>(cs[u].size); ++i) {
                    const auto v = static_cast<std::size_t>(cs[u].candidates[i]);
                    if (cs[u].marked[i] >= 0) {
                        bitmap::set_bit(temp_bitmap, v);
                    }
                }
                simd_bitmap::and_inplace(orbit_remaining_.row_mut(oid), temp_bitmap.data(), bitmap_words_);

                bitmap::clear(temp_bitmap);
                for (std::size_t i = 0; i < static_cast<std::size_t>(cs[u].size); ++i) {
                    const auto v = static_cast<std::size_t>(cs[u].candidates[i]);
                    if (cs[u].marked[i] == 1) {
                        bitmap::set_bit(temp_bitmap, v);
                    }
                }
                simd_bitmap::and_inplace(orbit_confirmed_.row_mut(oid), temp_bitmap.data(), bitmap_words_);
            }
        }

        orbit_remaining_count_[oid] = orbit_remaining_.popcount(oid);
        orbit_confirmed_count_[oid] = orbit_confirmed_.popcount(oid);
    }

    for (std::size_t oid = 0; oid < num_orbits; ++oid) {
        const auto start = orbit_members_offset_[oid];
        const auto size = orbit_size[oid];
        for (std::size_t idx = 0; idx < size; ++idx) {
            const auto u = orbit_members_data_[start + idx];
            if (u >= cs.size()) {
                continue;
            }
            for (std::size_t i = 0; i < static_cast<std::size_t>(cs[u].size); ++i) {
                const auto v = static_cast<std::size_t>(cs[u].candidates[i]);
                if (cs[u].marked[i] >= 0 && !orbit_remaining_.test_bit(oid, v)) {
                    cs[u].marked[i] = -1;
                }
            }
        }
    }

    orbit_lb_.assign(num_orbits, 0);
    orbit_ub_.assign(num_orbits, 0);
    for (std::size_t oid = 0; oid < num_orbits; ++oid) {
        const auto size = orbit_size[oid];
        if (size > 0) {
            orbit_lb_[oid] = static_cast<std::int32_t>(orbit_confirmed_count_[oid] / size);
            orbit_ub_[oid] = static_cast<std::int32_t>(orbit_remaining_count_[oid] / size);
        }
    }

    class_lb_.assign(num_labels, 0);
    class_ub_.assign(num_labels, std::numeric_limits<std::int32_t>::max());
    class_dirty_.assign(num_labels, 1);

    temp_union_.assign(bitmap_words_, 0ULL);
    temp_union2_.assign(bitmap_words_, 0ULL);
    total_verified = 0;
    total_success = 0;
    cs_cidx_.assign(static_cast<std::size_t>(q.n_vertex), 0);
    sync_buffer_.clear();
    class_subset_denom_.assign(num_labels, {});
    orbit_intersect_.assign(num_labels, {});
    orbit_index_in_label_.assign(num_orbits, 0);
    temp_intersect_.assign(bitmap_words_, 0ULL);
    orbit_has_unverified_.assign(num_orbits, 0);
    orbit_unverified_count_.assign(num_orbits, 0);
    vertex_unverified_count_.assign(static_cast<std::size_t>(q.n_vertex), 0);

    for (std::size_t label_idx = 0; label_idx < num_labels; ++label_idx) {
        const auto& orbit_ids = label_to_orbits[label_idx];
        const auto m = orbit_ids.size();
        if (m > 0) {
            const auto num_subsets = (static_cast<std::size_t>(1) << m) - 1;
            class_subset_denom_[label_idx].assign(num_subsets, 0);
            for (std::size_t mask = 1; mask <= num_subsets; ++mask) {
                std::size_t denom = 0;
                for (std::size_t i = 0; i < m; ++i) {
                    if (((mask >> i) & 1U) != 0U) {
                        denom += orbit_size[orbit_ids[i]];
                    }
                }
                class_subset_denom_[label_idx][mask - 1] = denom;
            }
        }
    }

    for (std::size_t label_idx = 0; label_idx < num_labels; ++label_idx) {
        const auto& orbit_ids = label_to_orbits[label_idx];
        const auto m = orbit_ids.size();
        for (std::size_t i = 0; i < m; ++i) {
            orbit_index_in_label_[orbit_ids[i]] = i;
        }
        if (m <= 1) {
            orbit_intersect_[label_idx].clear();
            continue;
        }

        auto intersect = std::vector<std::uint8_t>(m * m, 1);
        for (std::size_t i = 0; i < m; ++i) {
            intersect[i * m + i] = 1;
            for (std::size_t j = i + 1; j < m; ++j) {
                std::copy_n(orbit_remaining_.row(orbit_ids[i]), bitmap_words_, temp_intersect_.data());
                simd_bitmap::and_inplace(temp_intersect_.data(), orbit_remaining_.row(orbit_ids[j]), bitmap_words_);
                const auto has_intersect =
                    static_cast<std::uint8_t>(simd_bitmap::popcount(temp_intersect_.data(), bitmap_words_) > 0 ? 1 : 0);
                intersect[i * m + j] = has_intersect;
                intersect[j * m + i] = has_intersect;
            }
        }
        orbit_intersect_[label_idx] = std::move(intersect);
    }

    for (std::size_t oid = 0; oid < num_orbits; ++oid) {
        const auto start = orbit_members_offset_[oid];
        const auto size = orbit_size[oid];
        std::size_t orbit_total = 0;
        for (std::size_t idx = 0; idx < size; ++idx) {
            const auto u = orbit_members_data_[start + idx];
            if (u >= cs.size()) {
                continue;
            }
            std::size_t vertex_count = 0;
            for (std::size_t i = 0; i < static_cast<std::size_t>(cs[u].size); ++i) {
                if (cs[u].marked[i] == 0) {
                    ++vertex_count;
                }
            }
            vertex_unverified_count_[u] = vertex_count;
            orbit_total += vertex_count;
        }
        orbit_unverified_count_[oid] = orbit_total;
        orbit_has_unverified_[oid] = orbit_total > 0;
    }

    update_global_bounds();
}

void MDSContext::mark_confirmed(std::size_t u, std::int32_t v) {
    if (u >= orbits.size()) {
        return;
    }
    const auto oid = orbits[u];
    const auto v_idx = static_cast<std::size_t>(v);
    if (orbit_confirmed_.test_bit(oid, v_idx)) {
        return;
    }

    orbit_confirmed_.set_bit(oid, v_idx);
    ++orbit_confirmed_count_[oid];
    if (orbit_size[oid] > 0) {
        orbit_lb_[oid] = static_cast<std::int32_t>(orbit_confirmed_count_[oid] / orbit_size[oid]);
    }

    const auto label = orbit_label[oid];
    if (label >= 0 && static_cast<std::size_t>(label) < class_dirty_.size()) {
        class_dirty_[static_cast<std::size_t>(label)] = 1;
    }
}

void MDSContext::update_orbit_intersect_for_failed(std::size_t label_idx, std::size_t oid, std::size_t v_idx) {
    (void)v_idx;
    const auto& orbit_ids = label_to_orbits[label_idx];
    const auto m = orbit_ids.size();
    if (m <= 1) {
        return;
    }

    if (oid >= orbit_index_in_label_.size()) {
        return;
    }
    const auto oid_idx = orbit_index_in_label_[oid];
    if (oid_idx >= m) {
        return;
    }
    auto& intersect = orbit_intersect_[label_idx];
    if (intersect.empty()) {
        return;
    }

    for (std::size_t j = 0; j < m; ++j) {
        if (j == oid_idx) {
            continue;
        }
        const auto ij = oid_idx * m + j;
        if (intersect[ij] != 1) {
            continue;
        }
        std::copy_n(orbit_remaining_.row(oid), bitmap_words_, temp_intersect_.data());
        simd_bitmap::and_inplace(temp_intersect_.data(), orbit_remaining_.row(orbit_ids[j]), bitmap_words_);
        const auto still_intersect = simd_bitmap::popcount(temp_intersect_.data(), bitmap_words_) > 0;
        if (!still_intersect) {
            intersect[ij] = 0;
            intersect[j * m + oid_idx] = 0;
            class_dirty_[label_idx] = 1;
        }
    }
}

void MDSContext::mark_failed(std::size_t u, std::int32_t v) {
    if (u >= orbits.size()) {
        return;
    }
    const auto oid = orbits[u];
    const auto v_idx = static_cast<std::size_t>(v);
    if (!orbit_remaining_.test_bit(oid, v_idx)) {
        return;
    }

    orbit_remaining_.clear_bit(oid, v_idx);
    orbit_remaining_count_[oid] = orbit_remaining_count_[oid] == 0 ? 0 : orbit_remaining_count_[oid] - 1;
    if (orbit_size[oid] > 0) {
        orbit_ub_[oid] = static_cast<std::int32_t>(orbit_remaining_count_[oid] / orbit_size[oid]);
    }

    const auto label = orbit_label[oid];
    if (label >= 0 && static_cast<std::size_t>(label) < class_dirty_.size()) {
        class_dirty_[static_cast<std::size_t>(label)] = 1;
        update_orbit_intersect_for_failed(static_cast<std::size_t>(label), oid, v_idx);
    }
}

void MDSContext::sync_confirmed_to_cs(std::size_t u, std::int32_t v, std::vector<CandidateSpace>& cs) const {
    if (u >= orbits.size()) {
        return;
    }
    const auto oid = orbits[u];
    const auto start = orbit_members_offset_[oid];
    const auto size = orbit_size[oid];
    for (std::size_t idx = 0; idx < size; ++idx) {
        const auto uj = orbit_members_data_[start + idx];
        if (uj != u && uj < cs.size() && static_cast<std::size_t>(v) < cs[uj].vertex_to_pos.size()) {
            const auto pos = cs[uj].get_vertex_pos(static_cast<std::size_t>(v));
            if (pos >= 0 && static_cast<std::size_t>(pos) < static_cast<std::size_t>(cs[uj].size) &&
                cs[uj].marked[static_cast<std::size_t>(pos)] == 0) {
                cs[uj].marked[static_cast<std::size_t>(pos)] = 1;
            }
        }
    }
}

void MDSContext::sync_failed_to_cs(std::size_t u, std::int32_t v, std::vector<CandidateSpace>& cs) {
    sync_buffer_.clear();
    if (u >= orbits.size()) {
        return;
    }
    const auto oid = orbits[u];
    mark_failed(u, v);
    const auto start = orbit_members_offset_[oid];
    const auto size = orbit_size[oid];
    for (std::size_t idx = 0; idx < size; ++idx) {
        const auto uj = orbit_members_data_[start + idx];
        if (uj != u && uj < cs.size() && static_cast<std::size_t>(v) < cs[uj].vertex_to_pos.size()) {
            const auto pos = cs[uj].get_vertex_pos(static_cast<std::size_t>(v));
            if (pos >= 0 && static_cast<std::size_t>(pos) < static_cast<std::size_t>(cs[uj].size) &&
                cs[uj].marked[static_cast<std::size_t>(pos)] == 0) {
                cs[uj].marked[static_cast<std::size_t>(pos)] = -1;
                sync_buffer_.emplace_back(uj, v);
            }
        }
    }
}

const std::vector<std::pair<std::size_t, std::int32_t>>& MDSContext::get_synced_failed() const {
    return sync_buffer_;
}

void MDSContext::batch_sync_confirmed(const std::vector<CandidateSpace>& cs) {
    for (std::size_t u = 0; u < cs.size(); ++u) {
        if (u >= orbits.size()) {
            continue;
        }
        const auto oid = orbits[u];
        for (std::size_t i = 0; i < static_cast<std::size_t>(cs[u].size); ++i) {
            if (cs[u].marked[i] == 1) {
                const auto v = static_cast<std::size_t>(cs[u].candidates[i]);
                if (!orbit_confirmed_.test_bit(oid, v)) {
                    orbit_confirmed_.set_bit(oid, v);
                    ++orbit_confirmed_count_[oid];
                    const auto label = orbit_label[oid];
                    if (label >= 0 && static_cast<std::size_t>(label) < class_dirty_.size()) {
                        class_dirty_[static_cast<std::size_t>(label)] = 1;
                    }
                }
            }
        }
        if (orbit_size[oid] > 0) {
            orbit_lb_[oid] = static_cast<std::int32_t>(orbit_confirmed_count_[oid] / orbit_size[oid]);
        }
    }
}

std::vector<std::vector<std::size_t>> MDSContext::build_connected_components(
    std::size_t label_idx,
    const std::vector<std::size_t>& orbit_ids) const {
    const auto m = orbit_ids.size();
    if (m <= 1) {
        return m == 1 ? std::vector<std::vector<std::size_t>>{{0}} : std::vector<std::vector<std::size_t>>{};
    }

    const auto& intersect = orbit_intersect_[label_idx];
    if (intersect.empty()) {
        std::vector<std::size_t> full(m);
        for (std::size_t i = 0; i < m; ++i) {
            full[i] = i;
        }
        return {std::move(full)};
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
                if (!visited[v]) {
                    const auto has_edge = intersect[u * m + v] != 0;
                    if (has_edge) {
                        visited[v] = 1;
                        stack.push_back(v);
                    }
                }
            }
        }
        if (!component.empty()) {
            components.push_back(std::move(component));
        }
    }
    return components;
}

std::pair<std::int32_t, std::int32_t> MDSContext::compute_component_bounds(
    const std::vector<std::size_t>& component,
    const std::vector<std::size_t>& orbit_ids,
    std::size_t label_idx) {
    (void)label_idx;
    const auto m = component.size();
    if (m == 0) {
        return {0, 0};
    }
    if (m == 1) {
        const auto oid = orbit_ids[component[0]];
        return {orbit_lb_[oid], orbit_ub_[oid]};
    }

    std::int32_t best_lb = std::numeric_limits<std::int32_t>::max();
    std::int32_t best_ub = std::numeric_limits<std::int32_t>::max();
    const auto num_subsets = (static_cast<std::size_t>(1) << m) - 1;
    std::vector<std::size_t> denoms(num_subsets, 0);
    for (std::size_t mask = 1; mask <= num_subsets; ++mask) {
        std::size_t denom = 0;
        for (std::size_t i = 0; i < m; ++i) {
            if (((mask >> i) & 1U) != 0U) {
                denom += orbit_size[orbit_ids[component[i]]];
            }
        }
        denoms[mask - 1] = denom;
    }

    for (std::uint32_t mask = 1; mask < (1U << m); ++mask) {
        const auto total_size = denoms[mask - 1];
        if (total_size == 0) {
            continue;
        }
        if (static_cast<std::int32_t>(total_size) >= best_ub && static_cast<std::int32_t>(total_size) >= best_lb) {
            continue;
        }

        std::fill(temp_union_.begin(), temp_union_.end(), 0ULL);
        for (std::size_t i = 0; i < m; ++i) {
            if (((mask >> i) & 1U) != 0U) {
                simd_bitmap::or_inplace(temp_union_.data(), orbit_remaining_.row(orbit_ids[component[i]]), bitmap_words_);
            }
        }
        const auto ub = static_cast<std::int32_t>(simd_bitmap::popcount(temp_union_.data(), bitmap_words_) / total_size);
        best_ub = std::min(best_ub, ub);
        if (best_ub <= tau) {
            best_lb = best_lb == std::numeric_limits<std::int32_t>::max() ? 0 : best_lb;
            break;
        }

        std::fill(temp_union2_.begin(), temp_union2_.end(), 0ULL);
        for (std::size_t i = 0; i < m; ++i) {
            if (((mask >> i) & 1U) != 0U) {
                simd_bitmap::or_inplace(temp_union2_.data(), orbit_confirmed_.row(orbit_ids[component[i]]), bitmap_words_);
            }
        }
        const auto lb =
            static_cast<std::int32_t>(simd_bitmap::popcount(temp_union2_.data(), bitmap_words_) / total_size);
        best_lb = std::min(best_lb, lb);
    }

    if (best_lb == std::numeric_limits<std::int32_t>::max()) {
        best_lb = 0;
    }
    if (best_ub == std::numeric_limits<std::int32_t>::max()) {
        best_ub = 0;
    }
    return {best_lb, best_ub};
}

std::pair<std::int32_t, std::int32_t> MDSContext::compute_class_bounds(std::size_t label_idx) {
    if (!class_dirty_[label_idx]) {
        return {class_lb_[label_idx], class_ub_[label_idx]};
    }

    const auto& orbit_ids = label_to_orbits[label_idx];
    const auto m = orbit_ids.size();
    if (m == 0) {
        class_lb_[label_idx] = 0;
        class_ub_[label_idx] = 0;
        class_dirty_[label_idx] = 0;
        return {0, 0};
    }
    if (m == 1) {
        const auto oid = orbit_ids[0];
        class_lb_[label_idx] = orbit_lb_[oid];
        class_ub_[label_idx] = orbit_ub_[oid];
        class_dirty_[label_idx] = 0;
        return {class_lb_[label_idx], class_ub_[label_idx]};
    }

    const auto components = build_connected_components(label_idx, orbit_ids);
    std::int32_t best_lb = std::numeric_limits<std::int32_t>::max();
    std::int32_t best_ub = std::numeric_limits<std::int32_t>::max();
    for (const auto& component : components) {
        const auto [comp_lb, comp_ub] = compute_component_bounds(component, orbit_ids, label_idx);
        best_lb = std::min(best_lb, comp_lb);
        best_ub = std::min(best_ub, comp_ub);
        if (best_ub <= tau) {
            class_lb_[label_idx] = best_lb == std::numeric_limits<std::int32_t>::max() ? 0 : best_lb;
            class_ub_[label_idx] = best_ub;
            class_dirty_[label_idx] = 0;
            return {class_lb_[label_idx], best_ub};
        }
    }

    if (best_lb == std::numeric_limits<std::int32_t>::max()) {
        best_lb = 0;
    }
    if (best_ub == std::numeric_limits<std::int32_t>::max()) {
        best_ub = 0;
    }
    class_lb_[label_idx] = best_lb;
    class_ub_[label_idx] = best_ub;
    class_dirty_[label_idx] = 0;
    return {best_lb, best_ub};
}

void MDSContext::update_global_bounds() {
    std::int32_t lb = std::numeric_limits<std::int32_t>::max();
    std::int32_t ub = std::numeric_limits<std::int32_t>::max();
    for (const auto label_idx : active_labels) {
        const auto [clb, cub] = compute_class_bounds(label_idx);
        lb = std::min(lb, clb);
        ub = std::min(ub, cub);
    }
    global_lb = lb == std::numeric_limits<std::int32_t>::max() ? 0 : lb;
    global_ub = ub == std::numeric_limits<std::int32_t>::max() ? 0 : ub;
}

bool MDSContext::is_converged() const {
    return global_lb >= global_ub;
}

bool MDSContext::can_prune(std::int32_t tau_in) const {
    return global_ub <= tau_in;
}

std::int32_t MDSContext::get_mds() const {
    return global_lb;
}

double MDSContext::success_rate() const {
    if (total_verified == 0) {
        return 0.5;
    }
    return static_cast<double>(total_success) / static_cast<double>(total_verified);
}

std::size_t MDSContext::get_orbit_remaining_count(std::size_t u) const {
    if (u < orbits.size()) {
        return orbit_remaining_count_[orbits[u]];
    }
    return 0;
}

std::optional<std::pair<std::size_t, std::size_t>> MDSContext::select_for_raise_lb(
    const std::vector<CandidateSpace>& cs) {
    std::int32_t min_lb = std::numeric_limits<std::int32_t>::max();
    std::optional<std::size_t> min_label;
    for (const auto label_idx : active_labels) {
        if (class_lb_[label_idx] < min_lb) {
            min_lb = class_lb_[label_idx];
            min_label = label_idx;
        }
    }
    if (!min_label.has_value()) {
        return std::nullopt;
    }

    for (const auto oid : label_to_orbits[*min_label]) {
        if (!orbit_has_unverified_[oid]) {
            continue;
        }
        const auto member_start = orbit_members_offset_[oid];
        const auto member_size = orbit_size[oid];
        for (std::size_t i = 0; i < member_size; ++i) {
            const auto u = orbit_members_data_[member_start + i];
            if (u >= cs.size() || u >= cs_cidx_.size()) {
                continue;
            }
            if (vertex_unverified_count_[u] == 0) {
                continue;
            }
            const auto cidx_start = cs_cidx_[u];
            for (std::size_t ci = cidx_start; ci < static_cast<std::size_t>(cs[u].size); ++ci) {
                if (cs[u].marked[ci] == 0) {
                    return std::make_optional(std::make_pair(u, ci));
                }
            }
            vertex_unverified_count_[u] = 0;
        }
        orbit_has_unverified_[oid] = 0;
        orbit_unverified_count_[oid] = 0;
    }
    return std::nullopt;
}

std::optional<std::pair<std::size_t, std::size_t>> MDSContext::select_for_lower_ub(
    const std::vector<CandidateSpace>& cs) {
    std::int32_t min_ub = std::numeric_limits<std::int32_t>::max();
    std::optional<std::size_t> min_label;
    for (const auto label_idx : active_labels) {
        if (class_ub_[label_idx] < min_ub) {
            min_ub = class_ub_[label_idx];
            min_label = label_idx;
        }
    }
    if (!min_label.has_value()) {
        return std::nullopt;
    }

    for (const auto oid : label_to_orbits[*min_label]) {
        if (!orbit_has_unverified_[oid]) {
            continue;
        }
        const auto member_start = orbit_members_offset_[oid];
        const auto member_size = orbit_size[oid];
        for (std::size_t i = 0; i < member_size; ++i) {
            const auto u = orbit_members_data_[member_start + i];
            if (u >= cs.size() || u >= cs_cidx_.size()) {
                continue;
            }
            if (vertex_unverified_count_[u] == 0) {
                continue;
            }
            const auto cidx_start = cs_cidx_[u];
            for (std::size_t ci = cidx_start; ci < static_cast<std::size_t>(cs[u].size); ++ci) {
                if (cs[u].marked[ci] == 0) {
                    return std::make_optional(std::make_pair(u, ci));
                }
            }
            vertex_unverified_count_[u] = 0;
        }
        orbit_has_unverified_[oid] = 0;
        orbit_unverified_count_[oid] = 0;
    }
    return std::nullopt;
}

std::optional<std::pair<std::size_t, std::size_t>> MDSContext::select_next(const std::vector<CandidateSpace>& cs) {
    if (success_rate() > threshold) {
        auto target = select_for_raise_lb(cs);
        if (target.has_value()) {
            return target;
        }
        return select_for_lower_ub(cs);
    }
    auto target = select_for_lower_ub(cs);
    if (target.has_value()) {
        return target;
    }
    return select_for_raise_lb(cs);
}

void MDSContext::advance_cidx(std::size_t u, std::size_t ci) {
    if (u < cs_cidx_.size()) {
        cs_cidx_[u] = ci + 1;
    }
    if (u < vertex_unverified_count_.size() && vertex_unverified_count_[u] > 0) {
        --vertex_unverified_count_[u];
        if (u < orbits.size()) {
            const auto oid = orbits[u];
            if (oid < orbit_unverified_count_.size() && orbit_unverified_count_[oid] > 0) {
                --orbit_unverified_count_[oid];
                if (orbit_unverified_count_[oid] == 0) {
                    orbit_has_unverified_[oid] = 0;
                }
            }
        }
    }
}

}
