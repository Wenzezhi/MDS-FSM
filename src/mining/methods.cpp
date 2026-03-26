#include "mds/mining/core.hpp"

#include <algorithm>
#include <limits>
#include <optional>
#include <utility>

#include "mds/dfs_code.hpp"
#include "mds/graph.hpp"
#include "mds/minimality_cache.hpp"
#include "mds/mining/heap.hpp"
#include "mds/mining/mds_context.hpp"
#include "mds/mining/orbit.hpp"
#include "mds/profiler.hpp"
#include "mds/simd_nlf.hpp"

namespace mds {

Graph MiningContext::create_subgraph_edge(std::int32_t l1, std::int32_t el, std::int32_t l2) const {
    profile::ScopedPhase phase("graph_build");
    Graph g;
    g.init_vertex(2);
    g.label[0] = l1;
    g.label[1] = l2;
    g.init_edge(1);
    g.degree = {1, 1};
    g.dfs_code.emplace_back(0, 1, l1, el, l2);
    g.right_most_path = {1, 0};
    g.rmp_size = 2;
    g.nbr_offset = {0, 1, 2};
    g.core = {1, 1};
    g.max_degree = 1;
    g.nbr[0] = {1, el};
    g.nbr[1] = {0, el};
    g.init_nbr_to_pos();
    return g;
}

Graph MiningContext::create_subgraph(const Graph& p,
                                     bool fwd,
                                     std::int32_t new_lbl,
                                     std::int32_t el,
                                     std::int32_t fst,
                                     std::int32_t snd) const {
    profile::ScopedPhase phase("graph_build");
    const auto nv = p.n_vertex + (fwd ? 1 : 0);
    Graph g;
    g.init_vertex(nv);
    for (std::size_t i = 0; i < static_cast<std::size_t>(p.n_vertex); ++i) {
        g.label[i] = p.label[i];
    }
    if (fwd) {
        g.label[static_cast<std::size_t>(nv - 1)] = new_lbl;
    }

    g.init_edge(p.n_edge + 1);
    for (std::size_t i = 0; i < static_cast<std::size_t>(p.n_vertex); ++i) {
        g.degree[i] = p.degree[i];
    }
    ++g.degree[static_cast<std::size_t>(fst)];
    ++g.degree[static_cast<std::size_t>(snd)];

    g.dfs_code = p.dfs_code;
    if (fwd) {
        g.dfs_code.emplace_back(fst, snd, g.label[static_cast<std::size_t>(fst)], el, new_lbl);
    } else {
        g.dfs_code.emplace_back(
            fst, snd, g.label[static_cast<std::size_t>(fst)], el, g.label[static_cast<std::size_t>(snd)]);
    }

    g.rmp_size = build_right_most_path(g.dfs_code, g.right_most_path);
    std::int32_t pos = 0;
    for (std::size_t i = 0; i < static_cast<std::size_t>(g.n_vertex); ++i) {
        g.nbr_offset[i] = pos;
        pos += g.degree[i];
        g.core[i] = g.degree[i];
        g.max_degree = std::max(g.max_degree, g.degree[i]);
    }
    g.nbr_offset[static_cast<std::size_t>(g.n_vertex)] = pos;

    std::size_t cp = 0;
    for (std::size_t i = 0; i < static_cast<std::size_t>(p.n_vertex); ++i) {
        for (std::size_t j = static_cast<std::size_t>(p.nbr_offset[i]);
             j < static_cast<std::size_t>(p.nbr_offset[i + 1]);
             ++j) {
            g.nbr[cp++] = p.nbr[j];
        }
        if (i == static_cast<std::size_t>(fst)) {
            g.nbr[cp++] = {snd, el};
        } else if (i == static_cast<std::size_t>(snd)) {
            g.nbr[cp++] = {fst, el};
        }
    }
    if (fwd) {
        g.nbr[cp] = {fst, el};
    }
    g.init_nbr_to_pos();
    return g;
}

bool MiningContext::is_in_freq_edge(std::int32_t l1, std::int32_t el, std::int32_t l2) const {
    const auto start = static_cast<std::size_t>(freq_edge_offset[static_cast<std::size_t>(l1)]);
    const auto end = static_cast<std::size_t>(freq_edge_offset[static_cast<std::size_t>(l1 + 1)]);
    for (std::size_t i = start; i < end; ++i) {
        if (std::get<1>(freq_edge[i]) == el && std::get<2>(freq_edge[i]) == l2) {
            return std::get<0>(freq_edge[i]) > tau;
        }
    }
    return false;
}

bool MiningContext::is_min(const Graph& g) {
    profile::ScopedPhase phase("is_min");
    const auto pattern_hash = MinimalityCache::hash_dfs_code(g.dfs_code);
    return minimality_cache.check(pattern_hash, [&]() { return is_min_uncached(g); });
}

bool MiningContext::is_min_uncached(const Graph& g) {
    const auto n = static_cast<std::size_t>(g.n_vertex);
    std::int32_t min_l = std::numeric_limits<std::int32_t>::max();
    for (std::size_t i = 0; i < n; ++i) {
        min_l = std::min(min_l, g.label[i]);
    }
    for (std::size_t i = 0; i < n; ++i) {
        if (g.label[i] == min_l) {
            auto rn = temp_vec_pool.acquire();
            auto ri = temp_vec_pool.acquire();
            auto stk = temp_vec_pool.acquire();
            rn.assign(n, -1);
            ri.assign(n, -1);
            stk.clear();
            rn[i] = 0;
            ri[0] = static_cast<std::int32_t>(i);
            stk.push_back(static_cast<std::int32_t>(i));
            const auto is_minimal = build_cmp_dfs(g, rn, ri, 1, 0, stk);
            release_vec(temp_vec_pool, std::move(stk));
            release_vec(temp_vec_pool, std::move(ri));
            release_vec(temp_vec_pool, std::move(rn));
            if (!is_minimal) {
                return false;
            }
        }
    }
    return true;
}

bool MiningContext::build_cmp_dfs(const Graph& g,
                                  std::vector<std::int32_t>& rn,
                                  std::vector<std::int32_t>& ri,
                                  std::int32_t cn,
                                  std::size_t dl,
                                  std::vector<std::int32_t>& stk) const {
    if (stk.empty()) {
        return true;
    }

    const auto node = static_cast<std::size_t>(stk.back());
    std::pair<std::int32_t, std::int32_t> min_l{std::numeric_limits<std::int32_t>::max(),
                                                std::numeric_limits<std::int32_t>::max()};
    bool only_leaf = true;

    for (std::size_t i = static_cast<std::size_t>(g.nbr_offset[node]);
         i < static_cast<std::size_t>(g.nbr_offset[node + 1]);
         ++i) {
        const auto nb = static_cast<std::size_t>(g.nbr[i].first);
        if (rn[nb] == -1) {
            const std::pair<std::int32_t, std::int32_t> p{g.nbr[i].second, g.label[nb]};
            if (p < min_l) {
                min_l = p;
                only_leaf = g.degree[nb] == 1;
            } else if (p == min_l && g.degree[nb] != 1) {
                only_leaf = false;
            }
        }
    }

    if (min_l.first == std::numeric_limits<std::int32_t>::max()) {
        stk.pop_back();
        return build_cmp_dfs(g, rn, ri, cn, dl, stk);
    }

    if (only_leaf) {
        std::vector<std::size_t> same_label_leaves;
        for (std::size_t i = static_cast<std::size_t>(g.nbr_offset[node]);
             i < static_cast<std::size_t>(g.nbr_offset[node + 1]);
             ++i) {
            const auto nb = static_cast<std::size_t>(g.nbr[i].first);
            if (rn[nb] == -1) {
                const std::pair<std::int32_t, std::int32_t> p{g.nbr[i].second, g.label[nb]};
                if (p == min_l && g.degree[nb] == 1) {
                    same_label_leaves.push_back(nb);
                }
            }
        }

        auto temp_cn = cn;
        auto temp_dl = dl;
        for (const auto neighbor : same_label_leaves) {
            rn[neighbor] = temp_cn;
            ri[static_cast<std::size_t>(temp_cn)] = static_cast<std::int32_t>(neighbor);
            ++temp_cn;

            const DfsCode fc(
                rn[node], rn[neighbor], g.label[node], min_l.first, g.label[neighbor]);

            if (temp_dl < g.dfs_code.size()) {
                if (g.dfs_code[temp_dl] < fc) {
                    for (const auto nb : same_label_leaves) {
                        rn[nb] = -1;
                    }
                    return true;
                }
                if (fc < g.dfs_code[temp_dl]) {
                    for (const auto nb : same_label_leaves) {
                        rn[nb] = -1;
                    }
                    return false;
                }
            }
            ++temp_dl;
        }

        if (!build_cmp_dfs(g, rn, ri, temp_cn, temp_dl, stk)) {
            for (const auto nb : same_label_leaves) {
                rn[nb] = -1;
            }
            return false;
        }

        for (const auto nb : same_label_leaves) {
            rn[nb] = -1;
        }
        return true;
    }

    for (std::size_t i = static_cast<std::size_t>(g.nbr_offset[node]);
         i < static_cast<std::size_t>(g.nbr_offset[node + 1]);
         ++i) {
        const auto nb = static_cast<std::size_t>(g.nbr[i].first);
        const auto el = g.nbr[i].second;
        if (el == min_l.first && g.label[nb] == min_l.second && rn[nb] == -1) {
            rn[nb] = cn;
            ri[static_cast<std::size_t>(cn)] = static_cast<std::int32_t>(nb);
            const DfsCode fc(rn[node], rn[nb], g.label[node], el, g.label[nb]);

            auto dlen = dl;
            if (dlen < g.dfs_code.size()) {
                if (g.dfs_code[dlen] < fc) {
                    rn[nb] = -1;
                    ri[static_cast<std::size_t>(cn)] = -1;
                    return true;
                }
                if (fc < g.dfs_code[dlen]) {
                    rn[nb] = -1;
                    ri[static_cast<std::size_t>(cn)] = -1;
                    return false;
                }
            }
            ++dlen;

            for (std::size_t k = 0; k < static_cast<std::size_t>(cn); ++k) {
                const auto n = static_cast<std::size_t>(ri[k]);
                if (n != node) {
                    for (std::size_t j = static_cast<std::size_t>(g.nbr_offset[n]);
                         j < static_cast<std::size_t>(g.nbr_offset[n + 1]);
                         ++j) {
                        if (static_cast<std::size_t>(g.nbr[j].first) == nb) {
                            const DfsCode bc(rn[nb], rn[n], g.label[nb], g.nbr[j].second, g.label[n]);
                            if (dlen < g.dfs_code.size()) {
                                if (g.dfs_code[dlen] < bc) {
                                    rn[nb] = -1;
                                    ri[static_cast<std::size_t>(cn)] = -1;
                                    return true;
                                }
                                if (bc < g.dfs_code[dlen]) {
                                    rn[nb] = -1;
                                    ri[static_cast<std::size_t>(cn)] = -1;
                                    return false;
                                }
                            }
                            ++dlen;
                            break;
                        }
                    }
                }
            }

            stk.push_back(static_cast<std::int32_t>(nb));
            if (!build_cmp_dfs(g, rn, ri, cn + 1, dlen, stk)) {
                rn[nb] = -1;
                ri[static_cast<std::size_t>(cn)] = -1;
                if (!stk.empty()) {
                    stk.pop_back();
                }
                return false;
            }
            if (!stk.empty()) {
                stk.pop_back();
            }
            rn[nb] = -1;
            ri[static_cast<std::size_t>(cn)] = -1;
        }
    }

    return true;
}

std::vector<CandidateSpace> MiningContext::alloc_cs(const Graph& q) const {
    const auto nq = static_cast<std::size_t>(q.n_vertex);
    std::vector<CandidateSpace> cs(nq);
    for (std::size_t u = 0; u < nq; ++u) {
        cs[u].init(max_num_candidate, static_cast<std::size_t>(q.degree[u]), max_num_data_vertex);
        if (use_filter) {
            cs[u].init_nbr_set_cnt(max_num_candidate, max_degree);
        }
    }
    return cs;
}

bool MiningContext::filter_count(const Graph& q) const {
    if (q.n_vertex > data_graph.n_vertex || q.n_edge > data_graph.n_edge || q.max_degree > data_graph.max_degree) {
        return false;
    }
    for (const auto& [label, count] : q.label_frequency) {
        const auto it = data_graph.label_frequency.find(label);
        if (it == data_graph.label_frequency.end() || count > it->second) {
            return false;
        }
    }
    return true;
}

bool MiningContext::build_cs_edge(const Graph& q, std::vector<CandidateSpace>& cs) {
    const auto u1 = static_cast<std::size_t>(q.dfs_code[0].from);
    const auto u2 = static_cast<std::size_t>(q.dfs_code[0].to);
    const auto lu1 = q.dfs_code[0].from_label;
    const auto lu2 = q.dfs_code[0].to_label;
    const auto el = q.dfs_code[0].edge_label;
    const auto a = std::min(lu1, lu2);
    const auto b = std::max(lu1, lu2);
    const auto it = edge_list.find({a, b, el});
    if (it == edge_list.end()) {
        return false;
    }
    const auto eid = static_cast<std::size_t>(it->second);

    if (edge_img_bitmaps.has_value()) {
        const auto& bitmaps = *edge_img_bitmaps;
        const auto s1 = static_cast<std::size_t>(data_graph.vertex_offset[static_cast<std::size_t>(lu1)]);
        const auto e1 = static_cast<std::size_t>(data_graph.vertex_offset[static_cast<std::size_t>(lu1 + 1)]);
        for (std::size_t j = s1; j < e1; ++j) {
            const auto v = data_graph.vertex[j];
            if (bitmaps.contains(eid, 0, static_cast<std::size_t>(v))) {
                const auto p = static_cast<std::size_t>(cs[u1].size);
                cs[u1].candidates[p] = v;
                cs[u1].set_vertex_pos(static_cast<std::size_t>(v), static_cast<std::int32_t>(p));
                ++cs[u1].size;
            }
        }

        const auto s2 = static_cast<std::size_t>(data_graph.vertex_offset[static_cast<std::size_t>(lu2)]);
        const auto e2 = static_cast<std::size_t>(data_graph.vertex_offset[static_cast<std::size_t>(lu2 + 1)]);
        for (std::size_t j = s2; j < e2; ++j) {
            const auto v = data_graph.vertex[j];
            if (bitmaps.contains(eid, 1, static_cast<std::size_t>(v))) {
                const auto p = static_cast<std::size_t>(cs[u2].size);
                cs[u2].candidates[p] = v;
                cs[u2].set_vertex_pos(static_cast<std::size_t>(v), static_cast<std::int32_t>(p));
                ++cs[u2].size;
            }
        }
    } else {
        const auto s1 = static_cast<std::size_t>(data_graph.vertex_offset[static_cast<std::size_t>(lu1)]);
        const auto e1 = static_cast<std::size_t>(data_graph.vertex_offset[static_cast<std::size_t>(lu1 + 1)]);
        for (std::size_t j = s1; j < e1; ++j) {
            const auto v = data_graph.vertex[j];
            if (edge_img_sets[eid].first.contains(v)) {
                const auto p = static_cast<std::size_t>(cs[u1].size);
                cs[u1].candidates[p] = v;
                cs[u1].set_vertex_pos(static_cast<std::size_t>(v), static_cast<std::int32_t>(p));
                ++cs[u1].size;
            }
        }

        const auto s2 = static_cast<std::size_t>(data_graph.vertex_offset[static_cast<std::size_t>(lu2)]);
        const auto e2 = static_cast<std::size_t>(data_graph.vertex_offset[static_cast<std::size_t>(lu2 + 1)]);
        for (std::size_t j = s2; j < e2; ++j) {
            const auto v = data_graph.vertex[j];
            if (edge_img_sets[eid].second.contains(v)) {
                const auto p = static_cast<std::size_t>(cs[u2].size);
                cs[u2].candidates[p] = v;
                cs[u2].set_vertex_pos(static_cast<std::size_t>(v), static_cast<std::int32_t>(p));
                ++cs[u2].size;
            }
        }
    }

    return cs[u1].size > tau && cs[u2].size > tau;
}

std::int32_t MiningContext::compute_ub(const Graph& q, const std::vector<CandidateSpace>& cs) const {
    auto mni = static_cast<std::int32_t>(max_num_data_vertex);
    for (std::size_t u = 0; u < static_cast<std::size_t>(q.n_vertex); ++u) {
        mni = std::min(mni, cs[u].size);
    }
    return mni;
}

std::optional<std::pair<std::vector<CandidateSpace>, std::int32_t>> MiningContext::build_cs_from_edge(Graph& q) {
    profile::ScopedPhase phase("build_cs_edge");
    if (!filter_count(q)) {
        return std::nullopt;
    }
    q.compute_nlf(nlf_size, n_unique_label);
    q.compute_nbr_safety();
    auto cs = alloc_cs(q);

    if (q.n_vertex == 2) {
        if (!build_cs_edge(q, cs)) {
            return std::nullopt;
        }
    } else {
        const auto u1 = static_cast<std::size_t>(q.dfs_code[0].from);
        const auto u2 = static_cast<std::size_t>(q.dfs_code[0].to);
        const auto lu1 = q.dfs_code[0].from_label;
        const auto lu2 = q.dfs_code[0].to_label;
        const auto el0 = q.dfs_code[0].edge_label;
        const auto a = std::min(lu1, lu2);
        const auto b = std::max(lu1, lu2);
        const auto eid_it = edge_list.find({a, b, el0});
        if (eid_it == edge_list.end()) {
            return std::nullopt;
        }
        const auto eid = static_cast<std::size_t>(eid_it->second);

        if (edge_img_bitmaps.has_value()) {
            const auto& bitmaps = *edge_img_bitmaps;
            const auto s1 = static_cast<std::size_t>(data_graph.vertex_offset[static_cast<std::size_t>(lu1)]);
            const auto e1 = static_cast<std::size_t>(data_graph.vertex_offset[static_cast<std::size_t>(lu1 + 1)]);
            for (std::size_t j = s1; j < e1; ++j) {
                const auto v = data_graph.vertex[j];
                if (bitmaps.contains(eid, 0, static_cast<std::size_t>(v))) {
                    const auto p = static_cast<std::size_t>(cs[u1].size);
                    cs[u1].candidates[p] = v;
                    cs[u1].set_vertex_pos(static_cast<std::size_t>(v), static_cast<std::int32_t>(p));
                    ++cs[u1].size;
                }
            }

            const auto s2 = static_cast<std::size_t>(data_graph.vertex_offset[static_cast<std::size_t>(lu2)]);
            const auto e2 = static_cast<std::size_t>(data_graph.vertex_offset[static_cast<std::size_t>(lu2 + 1)]);
            for (std::size_t j = s2; j < e2; ++j) {
                const auto v = data_graph.vertex[j];
                if (bitmaps.contains(eid, 1, static_cast<std::size_t>(v))) {
                    const auto p = static_cast<std::size_t>(cs[u2].size);
                    cs[u2].candidates[p] = v;
                    cs[u2].set_vertex_pos(static_cast<std::size_t>(v), static_cast<std::int32_t>(p));
                    ++cs[u2].size;
                }
            }
        } else {
            const auto s1 = static_cast<std::size_t>(data_graph.vertex_offset[static_cast<std::size_t>(lu1)]);
            const auto e1 = static_cast<std::size_t>(data_graph.vertex_offset[static_cast<std::size_t>(lu1 + 1)]);
            for (std::size_t j = s1; j < e1; ++j) {
                const auto v = data_graph.vertex[j];
                if (edge_img_sets[eid].first.contains(v)) {
                    const auto p = static_cast<std::size_t>(cs[u1].size);
                    cs[u1].candidates[p] = v;
                    cs[u1].set_vertex_pos(static_cast<std::size_t>(v), static_cast<std::int32_t>(p));
                    ++cs[u1].size;
                }
            }

            const auto s2 = static_cast<std::size_t>(data_graph.vertex_offset[static_cast<std::size_t>(lu2)]);
            const auto e2 = static_cast<std::size_t>(data_graph.vertex_offset[static_cast<std::size_t>(lu2 + 1)]);
            for (std::size_t j = s2; j < e2; ++j) {
                const auto v = data_graph.vertex[j];
                if (edge_img_sets[eid].second.contains(v)) {
                    const auto p = static_cast<std::size_t>(cs[u2].size);
                    cs[u2].candidates[p] = v;
                    cs[u2].set_vertex_pos(static_cast<std::size_t>(v), static_cast<std::int32_t>(p));
                    ++cs[u2].size;
                }
            }
        }

        const auto cv = static_cast<std::size_t>(q.n_vertex - 1);
        const auto cl = q.label[cv];
        const auto nbr = static_cast<std::size_t>(q.dfs_code.back().from);

        to_clean_index = 0;
        std::uint8_t check_val = 0;
        const auto cs_nbr_size = static_cast<std::size_t>(cs[nbr].size);
        for (std::size_t y = 0; y < cs_nbr_size; ++y) {
            const auto pc = cs[nbr].candidates[y];
            const auto [st, cnt] = data_graph.lookup_label_nbr_offset(pc, cl);
            if (cnt == 0) {
                continue;
            }
            const auto end = static_cast<std::size_t>(st + cnt);
            for (std::size_t z = static_cast<std::size_t>(st); z < end; ++z) {
                const auto cid = static_cast<std::size_t>(data_graph.nbr[z].first);
                if (is_visited[cid] == check_val) {
                    ++is_visited[cid];
                    if (check_val == 0) {
                        array_to_clean[to_clean_index++] = static_cast<std::int32_t>(cid);
                    }
                }
            }
        }
        ++check_val;

        const auto q_degree_cv = q.degree[cv];
        const auto q_core_cv = q.core[cv];
        const auto q_max_nbr_degree_cv = q.max_nbr_degree[cv];
        SimdNlfComparator comparator(nlf_size);
        for (std::size_t ci = 0; ci < to_clean_index; ++ci) {
            const auto cid = static_cast<std::size_t>(array_to_clean[ci]);
            if (is_visited[cid] != check_val) {
                continue;
            }

            const auto d_degree = data_graph.degree[cid];
            const auto d_core = data_graph.core[cid];
            const auto d_max_nbr = data_graph.max_nbr_degree[cid];
            if (d_degree < q_degree_cv || d_core < q_core_cv || d_max_nbr < q_max_nbr_degree_cv) {
                continue;
            }

            const auto* q_nlf_slice = q.nlf.data() + cv * nlf_size;
            const auto* d_nlf_slice = data_graph.nlf.data() + cid * nlf_size;
            if (comparator.check(q_nlf_slice, d_nlf_slice)) {
                const auto p = static_cast<std::size_t>(cs[cv].size);
                cs[cv].candidates[p] = static_cast<std::int32_t>(cid);
                cs[cv].set_vertex_pos(cid, static_cast<std::int32_t>(p));
                ++cs[cv].size;
            }
        }

        while (to_clean_index > 0) {
            --to_clean_index;
            is_visited[static_cast<std::size_t>(array_to_clean[to_clean_index])] = 0;
        }

        if (cs[cv].size <= tau) {
            return std::nullopt;
        }
    }

    if (use_filter) {
        build_nbr_cnt(q, cs);
        if (!cs_node_filtering(q, cs)) {
            return std::nullopt;
        }
    }

    const auto ub = compute_ub(q, cs);
    if (ub <= tau) {
        return std::nullopt;
    }
    return std::make_optional(std::make_pair(std::move(cs), ub));
}

std::optional<std::vector<CandidateSpace>> MiningContext::copy_cs_fwd(const Graph& q,
                                                                      const Graph& p,
                                                                      const std::vector<CandidateSpace>& pcs) {
    auto cs = alloc_cs(q);

    for (std::size_t i = 0; i < static_cast<std::size_t>(q.n_vertex - 1); ++i) {
        const auto pcs_i_size = static_cast<std::size_t>(pcs[i].size);
        for (std::size_t j = 0; j < pcs_i_size; ++j) {
            if (pcs[i].marked[j] >= 0) {
                const auto v = pcs[i].candidates[j];
                const auto dst_idx = static_cast<std::size_t>(cs[i].size);
                cs[i].candidates[dst_idx] = v;
                cs[i].set_vertex_pos(static_cast<std::size_t>(v), static_cast<std::int32_t>(dst_idx));

                if (use_filter) {
                    const auto pdeg = static_cast<std::size_t>(p.degree[i]);
                    const auto nbr_len = std::min({pdeg, cs[i].nbr_cnt_row_len(), pcs[i].nbr_cnt_row_len()});
                    std::copy_n(pcs[i].nbr_cnt.row(j), nbr_len, cs[i].nbr_cnt.row_mut(dst_idx));

                    const auto dv = static_cast<std::size_t>(data_graph.degree[static_cast<std::size_t>(v)]);
                    const auto set_len = std::min({dv, cs[i].nbr_set_cnt_row_len(), pcs[i].nbr_set_cnt_row_len()});
                    std::copy_n(pcs[i].nbr_set_cnt.row(j), set_len, cs[i].nbr_set_cnt.row_mut(dst_idx));

                    const auto safety_len =
                        std::min({pdeg, cs[i].nbr_safety_row_len(), pcs[i].nbr_safety_row_len()});
                    std::copy_n(pcs[i].nbr_safety.row(j), safety_len, cs[i].nbr_safety.row_mut(dst_idx));
                }
                ++cs[i].size;
            }
        }
    }

    const auto cv = static_cast<std::size_t>(q.n_vertex - 1);
    const auto cl = q.label[cv];
    const auto nbr_u = static_cast<std::size_t>(q.dfs_code.back().from);

    to_clean_index = 0;
    std::uint8_t check_val = 0;
    const auto cs_nbr_u_size = static_cast<std::size_t>(cs[nbr_u].size);
    for (std::size_t y = 0; y < cs_nbr_u_size; ++y) {
        const auto pc = cs[nbr_u].candidates[y];
        const auto [st, cnt] = data_graph.lookup_label_nbr_offset(pc, cl);
        if (cnt == 0) {
            continue;
        }
        for (std::size_t z = static_cast<std::size_t>(st); z < static_cast<std::size_t>(st + cnt); ++z) {
            const auto cid = static_cast<std::size_t>(data_graph.nbr[z].first);
            if (is_visited[cid] == check_val) {
                ++is_visited[cid];
                if (check_val == 0) {
                    array_to_clean[to_clean_index++] = static_cast<std::int32_t>(cid);
                }
            }
        }
    }
    ++check_val;

    const auto q_degree_cv = q.degree[cv];
    const auto q_core_cv = q.core[cv];
    const auto q_max_nbr_degree_cv = q.max_nbr_degree[cv];
    SimdNlfComparator comparator(nlf_size);
    for (std::size_t ci = 0; ci < to_clean_index; ++ci) {
        const auto cid = static_cast<std::size_t>(array_to_clean[ci]);
        if (is_visited[cid] != check_val) {
            continue;
        }

        const auto d_degree = data_graph.degree[cid];
        const auto d_core = data_graph.core[cid];
        const auto d_max_nbr = data_graph.max_nbr_degree[cid];
        if (d_degree < q_degree_cv || d_core < q_core_cv || d_max_nbr < q_max_nbr_degree_cv) {
            continue;
        }

        const auto* q_nlf_slice = q.nlf.data() + cv * nlf_size;
        const auto* d_nlf_slice = data_graph.nlf.data() + cid * nlf_size;
        if (comparator.check(q_nlf_slice, d_nlf_slice)) {
            const auto p_idx = static_cast<std::size_t>(cs[cv].size);
            cs[cv].candidates[p_idx] = static_cast<std::int32_t>(cid);
            cs[cv].set_vertex_pos(cid, static_cast<std::int32_t>(p_idx));
            ++cs[cv].size;
        }
    }

    while (to_clean_index > 0) {
        --to_clean_index;
        is_visited[static_cast<std::size_t>(array_to_clean[to_clean_index])] = 0;
    }

    if (cs[cv].size <= tau) {
        return std::nullopt;
    }

    if (use_filter) {
        const auto u1 = static_cast<std::size_t>(q.dfs_code.back().from);
        const auto u2 = static_cast<std::size_t>(q.dfs_code.back().to);
        const auto lu1 = q.dfs_code.back().from_label;
        const auto lu2 = q.dfs_code.back().to_label;
        const auto qel = q.dfs_code.back().edge_label;
        const auto u2pos = q.nbr_to_pos_at(u1, u2);
        const auto u1pos = q.nbr_to_pos_at(u2, u1);

        if (u2pos >= 0 && u1pos >= 0) {
            const auto u2pos_u = static_cast<std::size_t>(u2pos);
            const auto u1pos_u = static_cast<std::size_t>(u1pos);

            for (std::size_t j = 0; j < static_cast<std::size_t>(cs[u1].size); ++j) {
                const auto v = cs[u1].candidates[j];
                const auto [st, cnt] = data_graph.lookup_label_nbr_offset(v, lu2);
                if (cnt != 0) {
                    for (std::size_t z = static_cast<std::size_t>(st); z < static_cast<std::size_t>(st + cnt); ++z) {
                        const auto vn = data_graph.nbr[z].first;
                        const auto el = data_graph.nbr[z].second;
                        const auto vnpos = cs[u2].get_vertex_pos(static_cast<std::size_t>(vn));
                        if (vnpos >= 0 && qel == el) {
                            const auto vnpos_u = static_cast<std::size_t>(vnpos);
                            if (vnpos_u < cs[u2].nbr_cnt_len() && u1pos_u < cs[u2].nbr_cnt_row_len()) {
                                ++(*cs[u2].nbr_cnt.get_mut(vnpos_u, u1pos_u));
                            }

                            const auto idx =
                                z - static_cast<std::size_t>(data_graph.nbr_offset[static_cast<std::size_t>(v)]);
                            if (j < cs[u1].nbr_set_cnt_len() && idx < cs[u1].nbr_set_cnt_row_len()) {
                                ++(*cs[u1].nbr_set_cnt.get_mut(j, idx));
                                if (cs[u1].nbr_set_cnt.get(j, idx) == 1) {
                                    const auto comp_idx = static_cast<std::size_t>(q.comp_label_idx_at(u1, u2));
                                    if (j < cs[u1].nbr_safety_len() && comp_idx < cs[u1].nbr_safety_row_len()) {
                                        ++(*cs[u1].nbr_safety.get_mut(j, comp_idx));
                                    }
                                }
                            }
                        }
                    }
                }

                const auto comp_idx = static_cast<std::size_t>(q.comp_label_idx_at(u1, u2));
                if (j < cs[u1].nbr_safety_len() && comp_idx < cs[u1].nbr_safety_row_len()) {
                    --(*cs[u1].nbr_safety.get_mut(j, comp_idx));
                }
            }

            for (std::size_t j = 0; j < static_cast<std::size_t>(cs[u2].size); ++j) {
                const auto v = cs[u2].candidates[j];
                const auto [st, cnt] = data_graph.lookup_label_nbr_offset(v, lu1);
                if (cnt != 0) {
                    for (std::size_t z = static_cast<std::size_t>(st); z < static_cast<std::size_t>(st + cnt); ++z) {
                        const auto vn = data_graph.nbr[z].first;
                        const auto el = data_graph.nbr[z].second;
                        const auto vnpos = cs[u1].get_vertex_pos(static_cast<std::size_t>(vn));
                        if (vnpos >= 0 && qel == el) {
                            const auto vnpos_u = static_cast<std::size_t>(vnpos);
                            if (vnpos_u < cs[u1].nbr_cnt_len() && u2pos_u < cs[u1].nbr_cnt_row_len()) {
                                ++(*cs[u1].nbr_cnt.get_mut(vnpos_u, u2pos_u));
                            }

                            const auto idx =
                                z - static_cast<std::size_t>(data_graph.nbr_offset[static_cast<std::size_t>(v)]);
                            if (j < cs[u2].nbr_set_cnt_len() && idx < cs[u2].nbr_set_cnt_row_len()) {
                                ++(*cs[u2].nbr_set_cnt.get_mut(j, idx));
                                if (cs[u2].nbr_set_cnt.get(j, idx) == 1) {
                                    const auto comp_idx = static_cast<std::size_t>(q.comp_label_idx_at(u2, u1));
                                    if (j < cs[u2].nbr_safety_len() && comp_idx < cs[u2].nbr_safety_row_len()) {
                                        ++(*cs[u2].nbr_safety.get_mut(j, comp_idx));
                                    }
                                }
                            }
                        }
                    }
                }

                const auto comp_idx = static_cast<std::size_t>(q.comp_label_idx_at(u2, u1));
                if (j < cs[u2].nbr_safety_len() && comp_idx < cs[u2].nbr_safety_row_len()) {
                    --(*cs[u2].nbr_safety.get_mut(j, comp_idx));
                }
            }
        }
    }

    return cs;
}

std::optional<std::vector<CandidateSpace>> MiningContext::copy_cs_bwd(const Graph& q,
                                                                      const Graph& p,
                                                                      const std::vector<CandidateSpace>& pcs) {
    auto cs = alloc_cs(q);

    for (std::size_t i = 0; i < static_cast<std::size_t>(q.n_vertex); ++i) {
        const auto pcs_i_size = static_cast<std::size_t>(pcs[i].size);
        for (std::size_t j = 0; j < pcs_i_size; ++j) {
            if (pcs[i].marked[j] >= 0) {
                const auto v = pcs[i].candidates[j];
                const auto dst_idx = static_cast<std::size_t>(cs[i].size);
                cs[i].candidates[dst_idx] = v;
                cs[i].set_vertex_pos(static_cast<std::size_t>(v), static_cast<std::int32_t>(dst_idx));

                if (use_filter) {
                    const auto pdeg = static_cast<std::size_t>(p.degree[i]);
                    const auto nbr_len = std::min({pdeg, cs[i].nbr_cnt_row_len(), pcs[i].nbr_cnt_row_len()});
                    std::copy_n(pcs[i].nbr_cnt.row(j), nbr_len, cs[i].nbr_cnt.row_mut(dst_idx));

                    const auto dv = static_cast<std::size_t>(data_graph.degree[static_cast<std::size_t>(v)]);
                    const auto set_len = std::min({dv, cs[i].nbr_set_cnt_row_len(), pcs[i].nbr_set_cnt_row_len()});
                    std::copy_n(pcs[i].nbr_set_cnt.row(j), set_len, cs[i].nbr_set_cnt.row_mut(dst_idx));

                    const auto safety_len =
                        std::min({pdeg, cs[i].nbr_safety_row_len(), pcs[i].nbr_safety_row_len()});
                    std::copy_n(pcs[i].nbr_safety.row(j), safety_len, cs[i].nbr_safety.row_mut(dst_idx));
                }
                ++cs[i].size;
            }
        }
    }

    if (use_filter) {
        const auto u1 = static_cast<std::size_t>(q.dfs_code.back().from);
        const auto u2 = static_cast<std::size_t>(q.dfs_code.back().to);
        const auto lu1 = q.dfs_code.back().from_label;
        const auto lu2 = q.dfs_code.back().to_label;
        const auto qel = q.dfs_code.back().edge_label;
        const auto u2pos = q.nbr_to_pos_at(u1, u2);
        const auto u1pos = q.nbr_to_pos_at(u2, u1);

        if (u2pos >= 0 && u1pos >= 0) {
            const auto u2pos_u = static_cast<std::size_t>(u2pos);
            const auto u1pos_u = static_cast<std::size_t>(u1pos);

            for (std::size_t j = 0; j < static_cast<std::size_t>(cs[u1].size); ++j) {
                const auto v = cs[u1].candidates[j];
                const auto [st, cnt] = data_graph.lookup_label_nbr_offset(v, lu2);
                if (cnt != 0) {
                    for (std::size_t z = static_cast<std::size_t>(st); z < static_cast<std::size_t>(st + cnt); ++z) {
                        const auto vn = data_graph.nbr[z].first;
                        const auto el = data_graph.nbr[z].second;
                        const auto vnpos = cs[u2].get_vertex_pos(static_cast<std::size_t>(vn));
                        if (vnpos >= 0 && qel == el) {
                            const auto vnpos_u = static_cast<std::size_t>(vnpos);
                            if (vnpos_u < cs[u2].nbr_cnt_len() && u1pos_u < cs[u2].nbr_cnt_row_len()) {
                                ++(*cs[u2].nbr_cnt.get_mut(vnpos_u, u1pos_u));
                            }
                            const auto idx =
                                z - static_cast<std::size_t>(data_graph.nbr_offset[static_cast<std::size_t>(v)]);
                            if (j < cs[u1].nbr_set_cnt_len() && idx < cs[u1].nbr_set_cnt_row_len()) {
                                ++(*cs[u1].nbr_set_cnt.get_mut(j, idx));
                                if (cs[u1].nbr_set_cnt.get(j, idx) == 1) {
                                    const auto comp_idx = static_cast<std::size_t>(q.comp_label_idx_at(u1, u2));
                                    if (j < cs[u1].nbr_safety_len() && comp_idx < cs[u1].nbr_safety_row_len()) {
                                        ++(*cs[u1].nbr_safety.get_mut(j, comp_idx));
                                    }
                                }
                            }
                        }
                    }
                }
                const auto comp_idx = static_cast<std::size_t>(q.comp_label_idx_at(u1, u2));
                if (j < cs[u1].nbr_safety_len() && comp_idx < cs[u1].nbr_safety_row_len()) {
                    --(*cs[u1].nbr_safety.get_mut(j, comp_idx));
                }
            }

            for (std::size_t j = 0; j < static_cast<std::size_t>(cs[u2].size); ++j) {
                const auto v = cs[u2].candidates[j];
                const auto [st, cnt] = data_graph.lookup_label_nbr_offset(v, lu1);
                if (cnt != 0) {
                    for (std::size_t z = static_cast<std::size_t>(st); z < static_cast<std::size_t>(st + cnt); ++z) {
                        const auto vn = data_graph.nbr[z].first;
                        const auto el = data_graph.nbr[z].second;
                        const auto vnpos = cs[u1].get_vertex_pos(static_cast<std::size_t>(vn));
                        if (vnpos >= 0 && qel == el) {
                            const auto vnpos_u = static_cast<std::size_t>(vnpos);
                            if (vnpos_u < cs[u1].nbr_cnt_len() && u2pos_u < cs[u1].nbr_cnt_row_len()) {
                                ++(*cs[u1].nbr_cnt.get_mut(vnpos_u, u2pos_u));
                            }
                            const auto idx =
                                z - static_cast<std::size_t>(data_graph.nbr_offset[static_cast<std::size_t>(v)]);
                            if (j < cs[u2].nbr_set_cnt_len() && idx < cs[u2].nbr_set_cnt_row_len()) {
                                ++(*cs[u2].nbr_set_cnt.get_mut(j, idx));
                                if (cs[u2].nbr_set_cnt.get(j, idx) == 1) {
                                    const auto comp_idx = static_cast<std::size_t>(q.comp_label_idx_at(u2, u1));
                                    if (j < cs[u2].nbr_safety_len() && comp_idx < cs[u2].nbr_safety_row_len()) {
                                        ++(*cs[u2].nbr_safety.get_mut(j, comp_idx));
                                    }
                                }
                            }
                        }
                    }
                }
                const auto comp_idx = static_cast<std::size_t>(q.comp_label_idx_at(u2, u1));
                if (j < cs[u2].nbr_safety_len() && comp_idx < cs[u2].nbr_safety_row_len()) {
                    --(*cs[u2].nbr_safety.get_mut(j, comp_idx));
                }
            }
        }
    }

    return cs;
}

std::optional<std::pair<std::vector<CandidateSpace>, std::int32_t>> MiningContext::build_cs_fwd(
    Graph& q,
    const Graph& p,
    const std::vector<CandidateSpace>& pcs) {
    profile::ScopedPhase phase("build_cs_fwd_bwd");
    if (!filter_count(q)) {
        return std::nullopt;
    }
    q.compute_nlf(nlf_size, n_unique_label);
    q.compute_nbr_safety();
    auto cs = copy_cs_fwd(q, p, pcs);
    if (!cs.has_value()) {
        return std::nullopt;
    }
    if (use_filter && !cs_node_filtering(q, *cs)) {
        return std::nullopt;
    }
    const auto ub = compute_ub(q, *cs);
    if (ub <= tau) {
        return std::nullopt;
    }
    return std::make_optional(std::make_pair(std::move(*cs), ub));
}

std::optional<std::pair<std::vector<CandidateSpace>, std::int32_t>> MiningContext::build_cs_bwd(
    Graph& q,
    const Graph& p,
    const std::vector<CandidateSpace>& pcs) {
    profile::ScopedPhase phase("build_cs_fwd_bwd");
    if (!filter_count(q)) {
        return std::nullopt;
    }
    q.compute_nlf(nlf_size, n_unique_label);
    q.compute_nbr_safety();
    auto cs = copy_cs_bwd(q, p, pcs);
    if (!cs.has_value()) {
        return std::nullopt;
    }
    if (use_filter && !cs_node_filtering(q, *cs)) {
        return std::nullopt;
    }
    const auto ub = compute_ub(q, *cs);
    if (ub <= tau) {
        return std::nullopt;
    }
    return std::make_optional(std::make_pair(std::move(*cs), ub));
}

std::string MiningContext::format_pattern(const Graph& g) const {
    std::string result = "[";
    for (std::size_t i = 0; i < g.dfs_code.size(); ++i) {
        const auto& code = g.dfs_code[i];
        result += "(" + std::to_string(code.from) + "," + std::to_string(code.to) + "," +
                  std::to_string(code.edge_label) + ")";
        if (i + 1 != g.dfs_code.size()) {
            result += ",";
        }
    }
    result += "]";
    return result;
}

const std::vector<std::size_t>& MiningContext::get_or_compute_orbits(const Graph& q) {
    profile::ScopedPhase phase("orbit_total");
    const auto pattern_hash = hash_pattern(q.dfs_code);
    auto [it, inserted] = orbit_cache.try_emplace(pattern_hash);
    if (inserted || it->second.size() != static_cast<std::size_t>(q.n_vertex)) {
        it->second = compute_orbits(q);
    }
    return it->second;
}

std::uint64_t MiningContext::extend_edge(CandidateHeap& gs, const Graph& p, std::int32_t parent_mds) {
    (void)parent_mds;
    std::uint64_t valid_count = 0;
    std::uint64_t generated_count = 0;

    for (std::size_t rn = 0; rn < static_cast<std::size_t>(p.rmp_size); ++rn) {
        const auto node = p.right_most_path[rn];
        const auto min_l = p.label[static_cast<std::size_t>(p.right_most_path[static_cast<std::size_t>(p.rmp_size - 1)])];
        const auto nl = p.label[static_cast<std::size_t>(node)];
        const auto start = static_cast<std::size_t>(freq_edge_offset[static_cast<std::size_t>(nl)]);
        const auto end = static_cast<std::size_t>(freq_edge_offset[static_cast<std::size_t>(nl + 1)]);
        for (std::size_t fi = start; fi < end; ++fi) {
            const auto [em, el, new_l] = freq_edge[fi];
            ++generated_count;
            if (new_l < min_l || em <= tau) {
                continue;
            }

            auto g = create_subgraph(p, true, new_l, el, node, p.n_vertex);
            if (!is_min(g)) {
                continue;
            }

            auto built = build_cs_from_edge(g);
            if (built.has_value() && built->second > tau) {
                const auto& orbits = get_or_compute_orbits(g);
                const auto mds_ub = compute_mds_ub(g, built->first, orbits, *this);
                if (mds_ub > tau) {
                    {
                        profile::ScopedPhase push_phase("queue_push_total");
                        gs.push(HeapElement::from_payload(
                            mds_ub, g.n_vertex, false, std::move(g), std::move(built->first)));
                    }
                    ++valid_count;
                }
            }
        }
    }

    candidates_generated += generated_count;
    return valid_count;
}

std::uint64_t MiningContext::extend(CandidateHeap& gs,
                                    const Graph& p,
                                    const std::vector<CandidateSpace>& pcs,
                                    std::int32_t mni) {
    (void)mni;
    std::uint64_t valid_count = 0;
    std::uint64_t generated_count = 0;

    for (std::size_t rn = 0; rn < static_cast<std::size_t>(p.rmp_size); ++rn) {
        const auto node = p.right_most_path[rn];
        const auto min_l = p.label[static_cast<std::size_t>(p.right_most_path[static_cast<std::size_t>(p.rmp_size - 1)])];
        const auto nl = p.label[static_cast<std::size_t>(node)];
        const auto start = static_cast<std::size_t>(freq_edge_offset[static_cast<std::size_t>(nl)]);
        const auto end = static_cast<std::size_t>(freq_edge_offset[static_cast<std::size_t>(nl + 1)]);
        for (std::size_t fi = start; fi < end; ++fi) {
            const auto [em, el, new_l] = freq_edge[fi];
            ++generated_count;
            if (new_l < min_l) {
                continue;
            }
            if (em <= tau) {
                break;
            }

            auto g = create_subgraph(p, true, new_l, el, node, p.n_vertex);
            if (!is_min(g)) {
                continue;
            }

            auto built = build_cs_fwd(g, p, pcs);
            if (built.has_value() && built->second > tau) {
                const auto& orbits = get_or_compute_orbits(g);
                const auto mds_ub = compute_mds_ub(g, built->first, orbits, *this);
                if (mds_ub > tau) {
                    {
                        profile::ScopedPhase push_phase("queue_push_total");
                        gs.push(HeapElement::from_payload(
                            mds_ub, g.n_vertex, false, std::move(g), std::move(built->first)));
                    }
                    ++valid_count;
                }
            }
        }
    }

    for (std::size_t rn = static_cast<std::size_t>(p.rmp_size); rn-- > 1;) {
        const auto fst = p.right_most_path[0];
        const auto snd = p.right_most_path[rn];
        bool exists = false;
        for (std::size_t i = static_cast<std::size_t>(p.nbr_offset[static_cast<std::size_t>(snd)]);
             i < static_cast<std::size_t>(p.nbr_offset[static_cast<std::size_t>(snd) + 1]);
             ++i) {
            if (fst == p.nbr[i].first) {
                exists = true;
                break;
            }
        }
        if (exists) {
            continue;
        }
        for (std::int32_t el = 0; el < static_cast<std::int32_t>(n_unique_edge_label); ++el) {
            ++generated_count;
            if (!is_in_freq_edge(
                    p.label[static_cast<std::size_t>(fst)], el, p.label[static_cast<std::size_t>(snd)])) {
                continue;
            }

            auto g = create_subgraph(p, false, 0, el, fst, snd);
            if (!is_min(g)) {
                continue;
            }

            auto built = build_cs_bwd(g, p, pcs);
            if (built.has_value() && built->second > tau) {
                const auto& orbits = get_or_compute_orbits(g);
                const auto mds_ub = compute_mds_ub(g, built->first, orbits, *this);
                if (mds_ub > tau) {
                    {
                        profile::ScopedPhase push_phase("queue_push_total");
                        gs.push(HeapElement::from_payload(
                            mds_ub, g.n_vertex, true, std::move(g), std::move(built->first)));
                    }
                    ++valid_count;
                }
            }
        }
    }

    candidates_generated += generated_count;
    return valid_count;
}

bool MiningContext::prepare_mni(const Graph& q, std::vector<CandidateSpace>& cs) {
    profile::ScopedPhase phase("prepare_mni");
    for (std::size_t u = 0; u < static_cast<std::size_t>(q.n_vertex); ++u) {
        for (std::size_t j = 0; j < static_cast<std::size_t>(cs[u].size); ++j) {
            const auto cand = cs[u].candidates[j];
            cs[u].set_vertex_pos(static_cast<std::size_t>(cand), static_cast<std::int32_t>(j));
        }
    }

    for (std::size_t u = static_cast<std::size_t>(q.n_vertex); u-- > 0;) {
        const auto cl = q.label[u];
        for (std::size_t j = static_cast<std::size_t>(q.nbr_offset[u]);
             j < static_cast<std::size_t>(q.nbr_offset[u + 1]);
             ++j) {
            const auto un = static_cast<std::size_t>(q.nbr[j].first);
            const auto qel = q.nbr[j].second;
            const auto upos = q.nbr_to_pos_at(un, u);
            if (upos < 0) {
                continue;
            }
            const auto upos_u = static_cast<std::size_t>(upos);
            if (upos_u >= cs[un].adjacent.num_rows()) {
                continue;
            }
            for (std::size_t vp = 0; vp < static_cast<std::size_t>(cs[u].size); ++vp) {
                cs[un].adjacent.at_mut(upos_u, vp).clear();
            }
            for (std::size_t vnp = 0; vnp < static_cast<std::size_t>(cs[un].size); ++vnp) {
                const auto vn = cs[un].candidates[vnp];
                const auto [st, cnt] = data_graph.lookup_label_nbr_offset(vn, cl);
                if (cnt == 0) {
                    continue;
                }
                for (std::size_t z = static_cast<std::size_t>(st); z < static_cast<std::size_t>(st + cnt); ++z) {
                    const auto v = data_graph.nbr[z].first;
                    const auto el = data_graph.nbr[z].second;
                    if (el != qel) {
                        continue;
                    }
                    const auto vp = cs[u].get_vertex_pos(static_cast<std::size_t>(v));
                    if (vp >= 0) {
                        const auto vp_u = static_cast<std::size_t>(vp);
                        cs[un].adjacent.at_mut(upos_u, vp_u).push_back(static_cast<std::int32_t>(vnp));
                    }
                }
            }
        }
    }
    return true;
}

std::int32_t MiningContext::compute_mni(const Graph& q, std::vector<CandidateSpace>& cs) {
    profile::ScopedPhase phase("compute_mni_total");
    auto bt_ctx_local = std::move(bt_ctx);
    if (!bt_ctx_local) {
        return -1;
    }

    const auto nq = static_cast<std::size_t>(q.n_vertex);
    const auto& orbits = get_or_compute_orbits(q);
    bt_ctx_local->set_orbits(orbits);

    MDSContext mds_ctx;
    {
        profile::ScopedPhase init_phase("mds_init");
        mds_ctx.init(q, cs, orbits, static_cast<std::size_t>(data_graph.n_vertex));
    }
    mds_ctx.tau = tau;

    std::int64_t total_cu = 0;
    for (std::size_t i = 0; i < nq; ++i) {
        total_cu += cs[i].size;
    }
    const auto avg_ci = nq == 0 ? 0.0 : static_cast<double>(total_cu) / static_cast<double>(nq);
    mds_ctx.threshold = avg_ci > 0.0 ? (avg_ci - static_cast<double>(tau)) / avg_ci : 0.5;

    std::vector<std::pair<std::size_t, std::int32_t>> newly_failed;
    newly_failed.reserve(64);

    while (!mds_ctx.is_converged()) {
        if (mds_ctx.can_prune(tau)) {
            bt_ctx = std::move(bt_ctx_local);
            return -1;
        }

        const auto target = [&]() {
            profile::ScopedPhase select_phase("mds_select_next");
            return mds_ctx.select_next(cs);
        }();
        if (!target.has_value()) {
            break;
        }

        const auto [curr_u, ci] = *target;
        if (cs[curr_u].marked[ci] != 0) {
            mds_ctx.advance_cidx(curr_u, ci);
            continue;
        }

        const auto v = cs[curr_u].candidates[ci];
        mds_ctx.advance_cidx(curr_u, ci);
        ++mds_ctx.total_verified;

        const auto success = [&]() {
            profile::ScopedPhase backtrack_phase("backtrack_total");
            return bt_ctx_local->backtrack_once(q, cs, curr_u, ci);
        }();
        if (success || cs[curr_u].marked[ci] == 1) {
            profile::ScopedPhase sync_phase("sync_confirmed_failed");
            ++mds_ctx.total_success;
            const auto& newly_confirmed = bt_ctx_local->get_newly_confirmed();
            for (const auto& [ui, vi] : newly_confirmed) {
                mds_ctx.mark_confirmed(ui, vi);
            }
            mds_ctx.sync_confirmed_to_cs(curr_u, v, cs);
        } else {
            profile::ScopedPhase sync_phase("sync_confirmed_failed");
            cs[curr_u].marked[ci] = -1;
            mds_ctx.sync_failed_to_cs(curr_u, v, cs);

            if (use_filter) {
                auto cs_c_local = std::move(cs_c_buffer);
                if (cs_c_local.size() < nq) {
                    cs_c_local.resize(nq, 0);
                }
                for (std::size_t i = 0; i < nq; ++i) {
                    cs_c_local[i] = static_cast<std::int32_t>(mds_ctx.get_orbit_remaining_count(i));
                }

                newly_failed.clear();
                if (!cs_node_filtering2(q, cs, cs_c_local, curr_u, v, newly_failed)) {
                    cs_c_buffer = std::move(cs_c_local);
                    bt_ctx = std::move(bt_ctx_local);
                    return -1;
                }

                const auto& synced = mds_ctx.get_synced_failed();
                for (const auto& [uj, vj] : synced) {
                    if (!cs_node_filtering2(q, cs, cs_c_local, uj, vj, newly_failed)) {
                        cs_c_buffer = std::move(cs_c_local);
                        bt_ctx = std::move(bt_ctx_local);
                        return -1;
                    }
                }

                for (const auto& [ui, vi] : newly_failed) {
                    mds_ctx.mark_failed(ui, vi);
                }
                cs_c_buffer = std::move(cs_c_local);
            }
        }

        {
            profile::ScopedPhase bounds_phase("update_global_bounds");
            mds_ctx.update_global_bounds();
        }
    }

    const auto mds = mds_ctx.get_mds();
    bt_ctx = std::move(bt_ctx_local);
    return mds;
}

}
