#include "mds/mining/core.hpp"
#include "mds/profiler.hpp"

#include <algorithm>

namespace mds {

void MiningContext::build_nbr_cnt(const Graph& q, std::vector<CandidateSpace>& cs) const {
    profile::ScopedPhase phase("filter_total");
    profile::ScopedPhase subphase("build_nbr_cnt_total");
    for (std::size_t u = 0; u < static_cast<std::size_t>(q.n_vertex); ++u) {
        auto& cs_u = cs[u];
        const auto cs_u_size = static_cast<std::size_t>(cs_u.size);
        for (std::size_t j = 0; j < cs_u_size; ++j) {
            const auto v = cs_u.candidates[j];
            cs_u.set_vertex_pos(static_cast<std::size_t>(v), static_cast<std::int32_t>(j));
            std::fill_n(cs_u.nbr_set_cnt.row_mut(j), static_cast<std::size_t>(data_graph.degree[static_cast<std::size_t>(v)]), std::int8_t{0});
        }
    }

    for (std::size_t u = 0; u < static_cast<std::size_t>(q.n_vertex); ++u) {
        auto& cs_u = cs[u];
        const auto nbr_start = static_cast<std::size_t>(q.nbr_offset[u]);
        const auto nbr_end = static_cast<std::size_t>(q.nbr_offset[u + 1]);
        for (std::size_t j = nbr_start; j < nbr_end; ++j) {
            const auto child = static_cast<std::size_t>(q.nbr[j].first);
            const auto qel = q.nbr[j].second;
            const auto nbr_idx = j - nbr_start;
            const auto q_label_child = q.nbr_label_data[j];
            const auto comp_idx = static_cast<std::size_t>(q.nbr_comp_idx_data[j]);

            for (std::size_t v_idx = 0; v_idx < static_cast<std::size_t>(cs_u.size); ++v_idx) {
                const auto v = cs_u.candidates[v_idx];
                const auto [st, cnt] = data_graph.lookup_label_nbr_offset(v, q_label_child);
                if (cnt == 0) {
                    continue;
                }
                auto* const nbr_cnt_row = cs_u.nbr_cnt.row_mut(v_idx);
                auto* const nbr_set_row = cs_u.nbr_set_cnt.row_mut(v_idx);
                auto* const nbr_safety_row = cs_u.nbr_safety.row_mut(v_idx);
                const auto v_nbr_offset = static_cast<std::size_t>(data_graph.nbr_offset[static_cast<std::size_t>(v)]);
                for (auto z = static_cast<std::size_t>(st); z < static_cast<std::size_t>(st + cnt); ++z) {
                    const auto vn = data_graph.nbr[z].first;
                    const auto el = data_graph.nbr[z].second;
                    const auto vn_pos = cs[child].get_vertex_pos(static_cast<std::size_t>(vn));
                    if (vn_pos == -1 || el != qel) {
                        continue;
                    }
                    ++nbr_cnt_row[nbr_idx];
                    const auto idx = z - v_nbr_offset;
                    auto& set_cnt = nbr_set_row[idx];
                    ++set_cnt;
                    if (set_cnt == 1) {
                        ++nbr_safety_row[comp_idx];
                    }
                }
                --nbr_safety_row[comp_idx];
            }
        }
    }
}

void MiningContext::remove_node(std::size_t u, std::int32_t v, std::vector<CandidateSpace>& cs, const Graph& q) const {
    (void)q;
    const auto vpos_value = cs[u].get_vertex_pos(static_cast<std::size_t>(v));
    if (vpos_value < 0) {
        return;
    }
    const auto vpos = static_cast<std::size_t>(vpos_value);
    const auto last_idx = static_cast<std::size_t>(cs[u].size - 1);
    if (vpos != last_idx) {
        const auto y = cs[u].candidates[last_idx];
        cs[u].candidates[vpos] = y;
        cs[u].nbr_cnt.copy_row(vpos, last_idx);
        cs[u].nbr_set_cnt.copy_row(vpos, last_idx);
        cs[u].nbr_safety.copy_row(vpos, last_idx);
        cs[u].marked[vpos] = cs[u].marked[last_idx];
        cs[u].set_vertex_pos(static_cast<std::size_t>(y), static_cast<std::int32_t>(vpos));
    }
    cs[u].clear_vertex_pos(static_cast<std::size_t>(v));
    --cs[u].size;
}

bool MiningContext::cs_node_filtering(const Graph& q, std::vector<CandidateSpace>& cs) const {
    profile::ScopedPhase phase("filter_total");
    profile::ScopedPhase subphase("node_filter_total");
    if (!use_filter) {
        return true;
    }

    auto& deadnodes = filter_deadnodes_buffer;
    deadnodes.clear();
    auto& cs_bc = filter_cs_bc_buffer;
    if (cs_bc.size() < static_cast<std::size_t>(q.n_vertex)) {
        cs_bc.resize(static_cast<std::size_t>(q.n_vertex), 0);
    }
    const auto u1 = static_cast<std::size_t>(q.dfs_code.back().from);
    const auto u2 = static_cast<std::size_t>(q.dfs_code.back().to);
    std::size_t deadnode_head = 0;

    for (std::size_t i = 0; i < static_cast<std::size_t>(q.n_vertex); ++i) {
        cs_bc[i] = cs[i].size;
    }

    auto& cs_u1 = cs[u1];
    for (std::int32_t j = 0; j < cs[u1].size; ++j) {
        const auto ju = static_cast<std::size_t>(j);
        const auto v = cs_u1.candidates[ju];
        const auto* const nbr_cnt_row = cs_u1.nbr_cnt.row(ju);
        const auto* const nbr_safety_row = cs_u1.nbr_safety.row(ju);
        bool removed = false;
        for (auto k = static_cast<std::size_t>(q.nbr_offset[u1]); k < static_cast<std::size_t>(q.nbr_offset[u1 + 1]); ++k) {
            const auto nbr_idx = k - static_cast<std::size_t>(q.nbr_offset[u1]);
            const auto comp_idx = static_cast<std::size_t>(q.nbr_comp_idx_data[k]);
            if (nbr_cnt_row[nbr_idx] == 0 || nbr_safety_row[comp_idx] < 0) {
                remove_node(u1, v, cs, q);
                --cs_bc[u1];
                if (cs_bc[u1] <= tau) {
                    return false;
                }
                deadnodes.emplace_back(u1, v);
                --j;
                removed = true;
                break;
            }
        }
        if (removed) {
            continue;
        }
    }

    auto& cs_u2 = cs[u2];
    for (std::int32_t j = 0; j < cs[u2].size; ++j) {
        const auto ju = static_cast<std::size_t>(j);
        const auto v = cs_u2.candidates[ju];
        const auto* const nbr_cnt_row = cs_u2.nbr_cnt.row(ju);
        const auto* const nbr_safety_row = cs_u2.nbr_safety.row(ju);
        bool removed = false;
        for (auto k = static_cast<std::size_t>(q.nbr_offset[u2]); k < static_cast<std::size_t>(q.nbr_offset[u2 + 1]); ++k) {
            const auto nbr_idx = k - static_cast<std::size_t>(q.nbr_offset[u2]);
            const auto comp_idx = static_cast<std::size_t>(q.nbr_comp_idx_data[k]);
            if (nbr_cnt_row[nbr_idx] == 0 || nbr_safety_row[comp_idx] < 0) {
                remove_node(u2, v, cs, q);
                --cs_bc[u2];
                if (cs_bc[u2] <= tau) {
                    return false;
                }
                deadnodes.emplace_back(u2, v);
                --j;
                removed = true;
                break;
            }
        }
        if (removed) {
            continue;
        }
    }

    while (deadnode_head < deadnodes.size()) {
        const auto [u, v] = deadnodes[deadnode_head++];
        for (auto j = static_cast<std::size_t>(q.nbr_offset[u]); j < static_cast<std::size_t>(q.nbr_offset[u + 1]); ++j) {
            const auto un = static_cast<std::size_t>(q.nbr[j].first);
            const auto qel = q.nbr[j].second;
            const auto [st, cnt] = data_graph.lookup_label_nbr_offset(v, q.nbr_label_data[j]);
            if (cnt == 0) {
                continue;
            }
            const auto unpos_value = q.nbr_reverse_pos_data[j];
            if (unpos_value < 0) {
                continue;
            }
            const auto unpos = static_cast<std::size_t>(unpos_value);
            const auto comp_idx = static_cast<std::size_t>(q.nbr_reverse_comp_idx_data[j]);
            for (auto z = static_cast<std::size_t>(st); z < static_cast<std::size_t>(st + cnt); ++z) {
                const auto vn = data_graph.nbr[z].first;
                const auto el = data_graph.nbr[z].second;
                if (el != qel) {
                    continue;
                }
                const auto vnpos_value = cs[un].get_vertex_pos(static_cast<std::size_t>(vn));
                if (vnpos_value == -1) {
                    continue;
                }
                const auto vnpos = static_cast<std::size_t>(vnpos_value);
                auto& cs_un = cs[un];
                auto* const nbr_cnt_row = cs_un.nbr_cnt.row_mut(vnpos);
                if (--nbr_cnt_row[unpos] == 0) {
                    remove_node(un, vn, cs, q);
                    --cs_bc[un];
                    if (cs_bc[un] <= tau) {
                        return false;
                    }
                    deadnodes.emplace_back(un, vn);
                } else {
                    const auto nbridx = static_cast<std::size_t>(inv_nbr[z]);
                    auto* const nbr_set_row = cs_un.nbr_set_cnt.row_mut(vnpos);
                    if (--nbr_set_row[nbridx] == 0) {
                        auto* const nbr_safety_row = cs_un.nbr_safety.row_mut(vnpos);
                        if (--nbr_safety_row[comp_idx] < 0) {
                            remove_node(un, vn, cs, q);
                            --cs_bc[un];
                            if (cs_bc[un] <= tau) {
                                return false;
                            }
                            deadnodes.emplace_back(un, vn);
                        }
                    }
                }
            }
        }
    }
    return true;
}

bool MiningContext::cs_node_filtering2(const Graph& q,
                                       std::vector<CandidateSpace>& cs,
                                       std::vector<std::int32_t>& cs_c,
                                       std::size_t iu,
                                       std::int32_t iv,
                                       std::vector<std::pair<std::size_t, std::int32_t>>& newly_failed) const {
    profile::ScopedPhase phase("filter_total");
    profile::ScopedPhase subphase("node_filter_total");
    if (!use_filter) {
        return true;
    }

    auto& deadnodes = filter_deadnodes_buffer;
    deadnodes.clear();
    deadnodes.emplace_back(iu, iv);
    std::size_t deadnode_head = 0;

    while (deadnode_head < deadnodes.size()) {
        const auto [u, v] = deadnodes[deadnode_head++];
        for (auto j = static_cast<std::size_t>(q.nbr_offset[u]); j < static_cast<std::size_t>(q.nbr_offset[u + 1]); ++j) {
            const auto un = static_cast<std::size_t>(q.nbr[j].first);
            const auto qel = q.nbr[j].second;
            const auto [st, cnt] = data_graph.lookup_label_nbr_offset(v, q.nbr_label_data[j]);
            if (cnt == 0) {
                continue;
            }
            const auto unpos_value = q.nbr_reverse_pos_data[j];
            if (unpos_value < 0) {
                continue;
            }
            const auto unpos = static_cast<std::size_t>(unpos_value);
            const auto comp_idx = static_cast<std::size_t>(q.nbr_reverse_comp_idx_data[j]);
            for (auto z = static_cast<std::size_t>(st); z < static_cast<std::size_t>(st + cnt); ++z) {
                const auto vn = data_graph.nbr[z].first;
                const auto el = data_graph.nbr[z].second;
                if (el != qel) {
                    continue;
                }
                const auto vnpos_value = cs[un].get_vertex_pos(static_cast<std::size_t>(vn));
                if (vnpos_value == -1) {
                    continue;
                }
                const auto vnpos = static_cast<std::size_t>(vnpos_value);
                auto& cs_un = cs[un];
                if (cs_un.marked[vnpos] < 0) {
                    continue;
                }
                auto* const nbr_cnt_row = cs_un.nbr_cnt.row_mut(vnpos);
                if (--nbr_cnt_row[unpos] == 0) {
                    cs_un.marked[vnpos] = -1;
                    --cs_c[un];
                    newly_failed.emplace_back(un, vn);
                    if (cs_c[un] <= tau) {
                        return false;
                    }
                    deadnodes.emplace_back(un, vn);
                } else {
                    const auto nbridx = static_cast<std::size_t>(inv_nbr[z]);
                    auto* const nbr_set_row = cs_un.nbr_set_cnt.row_mut(vnpos);
                    if (--nbr_set_row[nbridx] == 0) {
                        auto* const nbr_safety_row = cs_un.nbr_safety.row_mut(vnpos);
                        if (--nbr_safety_row[comp_idx] < 0) {
                            cs_un.marked[vnpos] = -1;
                            --cs_c[un];
                            newly_failed.emplace_back(un, vn);
                            if (cs_c[un] <= tau) {
                                return false;
                            }
                            deadnodes.emplace_back(un, vn);
                        }
                    }
                }
            }
        }
    }
    return true;
}

}
