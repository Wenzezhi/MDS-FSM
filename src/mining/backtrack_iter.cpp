#include "mds/mining/backtrack_iter.hpp"

#include <bit>

#include "mds/profiler.hpp"

namespace mds {

BacktrackContext::BacktrackContext(std::size_t max_data_vertex,
                                   std::size_t max_candidate,
                                   std::size_t max_degree)
    : element(),
      queue(),
      ancestors(kMaxQueryVertex, {0ULL}),
      dyanc(kMaxQueryVertex, kSmallQueryDegree, 1),
      failing_set_size(1),
      is_redundant(false),
      mark1_vertex(std::min<std::size_t>(std::max<std::size_t>(max_candidate, 1), 10000), 0),
      map_to(max_data_vertex, -1),
      iec_data(kMaxQueryVertex * kSmallQueryDegree * std::min<std::size_t>(std::max<std::size_t>(max_candidate, 1), 10000), 0),
      iec_size(kMaxQueryVertex * kSmallQueryDegree, 0),
      iec_capacity(std::min<std::size_t>(std::max<std::size_t>(max_candidate, 1), 10000)) {
    (void)max_degree;
    element.reserve(kMaxQueryVertex);
    for (std::size_t i = 0; i < kMaxQueryVertex; ++i) {
        element.emplace_back(failing_set_size);
    }
    curr_mapping.fill(-1);
    cand_pos.fill(-1);
    n_mapped_parent.fill(0);
    cs_m.fill(0);
    cs_c.fill(0);
    cs_cidx.fill(0);
}

std::size_t BacktrackContext::iec_idx(std::size_t u, std::size_t p, std::size_t i) const {
    return (u * kSmallQueryDegree + p) * iec_capacity + i;
}

std::size_t BacktrackContext::iec_size_idx(std::size_t u, std::size_t p) {
    return u * kSmallQueryDegree + p;
}

std::int32_t BacktrackContext::get_iec_size(std::size_t u, std::size_t p) const {
    return iec_size[iec_size_idx(u, p)];
}

void BacktrackContext::set_iec_size(std::size_t u, std::size_t p, std::int32_t value) {
    iec_size[iec_size_idx(u, p)] = value;
}

std::int32_t BacktrackContext::get_iec(std::size_t u, std::size_t p, std::size_t i) const {
    return iec_data[iec_idx(u, p, i)];
}

void BacktrackContext::set_iec(std::size_t u, std::size_t p, std::size_t i, std::int32_t value) {
    iec_data[iec_idx(u, p, i)] = value;
}

void BacktrackContext::set_orbits(const std::vector<std::size_t>& new_orbits) {
    orbits = new_orbits;
    std::size_t num_orbits = 0;
    for (const auto orbit : orbits) {
        num_orbits = std::max(num_orbits, orbit + 1);
    }

    orbit_members_offset_.assign(num_orbits, 0);
    std::vector<std::size_t> orbit_sizes(num_orbits, 0);
    for (std::size_t u = 0; u < orbits.size(); ++u) {
        ++orbit_sizes[orbits[u]];
    }

    std::size_t total_size = 0;
    for (std::size_t oid = 0; oid < num_orbits; ++oid) {
        orbit_members_offset_[oid] = total_size;
        total_size += orbit_sizes[oid];
    }

    orbit_members_data_.assign(total_size, 0);
    auto current_pos = orbit_members_offset_;
    for (std::size_t u = 0; u < orbits.size(); ++u) {
        orbit_members_data_[current_pos[orbits[u]]++] = u;
    }
}

void BacktrackContext::sync_mark_confirmed(std::size_t u, std::int32_t v, std::vector<CandidateSpace>& cs) {
    if (orbits.empty() || u >= orbits.size()) {
        return;
    }
    const auto orbit_id = orbits[u];
    for (const auto uj : orbit_members_span(orbit_id)) {
        if (uj != u && uj < cs.size()) {
            const auto vj_pos = cs[uj].get_vertex_pos(static_cast<std::size_t>(v));
            if (vj_pos >= 0) {
                const auto pos = static_cast<std::size_t>(vj_pos);
                if (pos < static_cast<std::size_t>(cs[uj].size) && cs[uj].marked[pos] == 0) {
                    cs[uj].marked[pos] = 1;
                    ++cs_m[uj];
                }
            }
        }
    }
}

std::vector<std::pair<std::size_t, std::int32_t>> BacktrackContext::sync_mark_failed(
    std::size_t u,
    std::int32_t v,
    std::vector<CandidateSpace>& cs) {
    std::vector<std::pair<std::size_t, std::int32_t>> synced;
    if (orbits.empty() || u >= orbits.size()) {
        return synced;
    }
    const auto orbit_id = orbits[u];
    for (const auto uj : orbit_members_span(orbit_id)) {
        if (uj != u && uj < cs.size()) {
            const auto vj_pos = cs[uj].get_vertex_pos(static_cast<std::size_t>(v));
            if (vj_pos >= 0) {
                const auto pos = static_cast<std::size_t>(vj_pos);
                if (pos < static_cast<std::size_t>(cs[uj].size) && cs[uj].marked[pos] == 0) {
                    cs[uj].marked[pos] = -1;
                    --cs_c[uj];
                    synced.emplace_back(uj, v);
                }
            }
        }
    }
    return synced;
}

const std::vector<std::pair<std::size_t, std::int32_t>>& BacktrackContext::get_newly_confirmed() const {
    return newly_confirmed_;
}

void BacktrackContext::clear_newly_confirmed() {
    newly_confirmed_.clear();
}

std::span<const std::size_t> BacktrackContext::orbit_members_span(std::size_t oid) const {
    if (oid >= orbit_members_offset_.size()) {
        return {};
    }
    const auto start = orbit_members_offset_[oid];
    const auto end = (oid + 1 < orbit_members_offset_.size()) ? orbit_members_offset_[oid + 1] : orbit_members_data_.size();
    return {orbit_members_data_.data() + static_cast<std::ptrdiff_t>(start), end - start};
}

std::int64_t BacktrackContext::get_extendable_candidates(std::size_t u,
                                                         std::size_t cand_parent_idx,
                                                         std::size_t u_ngb_idx,
                                                         std::size_t num_parents_mapped,
                                                         const std::vector<CandidateSpace>& cs) {
    const auto iec_pos = num_parents_mapped - 1;
    set_iec_size(u, iec_pos, 0);

    if (num_parents_mapped == 1) {
        if (u_ngb_idx < cs[u].adjacent.num_rows() && cand_parent_idx < cs[u].adjacent.row_capacity()) {
            for (const auto cand_idx : cs[u].adjacent.at(u_ngb_idx, cand_parent_idx)) {
                const auto ci = static_cast<std::size_t>(cand_idx);
                if (ci < static_cast<std::size_t>(cs[u].size) && cs[u].marked[ci] >= 0) {
                    const auto sz = static_cast<std::size_t>(get_iec_size(u, 0));
                    if (sz < iec_capacity) {
                        set_iec(u, 0, sz, cand_idx);
                        set_iec_size(u, 0, static_cast<std::int32_t>(sz + 1));
                    }
                }
            }
        }
    } else {
        const auto prev_sz = static_cast<std::size_t>(get_iec_size(u, iec_pos - 1));
        if (u_ngb_idx < cs[u].adjacent.num_rows() && cand_parent_idx < cs[u].adjacent.row_capacity()) {
            const auto& adj = cs[u].adjacent.at(u_ngb_idx, cand_parent_idx);
            std::size_t j = 0;
            std::size_t k = 0;
            while (j < prev_sz && k < adj.size()) {
                const auto cand_idx = get_iec(u, iec_pos - 1, j);
                const auto new_cand_idx = adj[k];
                if (cand_idx == new_cand_idx) {
                    const auto sz = static_cast<std::size_t>(get_iec_size(u, iec_pos));
                    if (sz < iec_capacity) {
                        set_iec(u, iec_pos, sz, cand_idx);
                        set_iec_size(u, iec_pos, static_cast<std::int32_t>(sz + 1));
                    }
                    ++j;
                    ++k;
                } else if (cand_idx < new_cand_idx) {
                    ++j;
                } else {
                    ++k;
                }
            }
        }
    }
    return get_iec_size(u, iec_pos);
}

void BacktrackContext::update_ancestors(std::size_t u_curr, const Graph& q) {
    for (auto i = static_cast<std::size_t>(q.nbr_offset[u_curr]); i < static_cast<std::size_t>(q.nbr_offset[u_curr + 1]); ++i) {
        const auto uc = static_cast<std::size_t>(q.nbr[i].first);
        if (!is_mapped[uc]) {
            ++n_mapped_parent[uc];
            const auto u_num_par = static_cast<std::size_t>(n_mapped_parent[uc]);
            const auto up_num_par = static_cast<std::size_t>(n_mapped_parent[u_curr]);
            dyanc.add_ancestor(uc, u_curr, u_num_par, up_num_par);
        }
    }
}

void BacktrackContext::update_extendable_root(const Graph& q,
                                              const std::vector<CandidateSpace>& cs,
                                              std::size_t ri,
                                              std::size_t root_u) {
    for (auto j = static_cast<std::size_t>(q.nbr_offset[root_u]); j < static_cast<std::size_t>(q.nbr_offset[root_u + 1]); ++j) {
        const auto rc = static_cast<std::size_t>(q.nbr[j].first);
        const auto u_ngb_idx = q.nbr_to_pos_at(rc, root_u);
        if (u_ngb_idx < 0) {
            continue;
        }
        const auto np = static_cast<std::size_t>(n_mapped_parent[rc]);
        (void)get_extendable_candidates(rc, ri, static_cast<std::size_t>(u_ngb_idx), np, cs);
        if (n_mapped_parent[rc] == 1) {
            queue.insert_to_queue(static_cast<std::int32_t>(rc));
        }
    }
}

bool BacktrackContext::update_extendable_u(const Graph& q,
                                           const std::vector<CandidateSpace>& cs,
                                           std::size_t u_curr,
                                           std::size_t depth,
                                           std::size_t u_pos) {
    for (auto j = static_cast<std::size_t>(q.nbr_offset[u_curr]); j < static_cast<std::size_t>(q.nbr_offset[u_curr + 1]); ++j) {
        const auto uc = static_cast<std::size_t>(q.nbr[j].first);
        if (is_mapped[uc]) {
            continue;
        }
        const auto u_ngb_idx = q.nbr_to_pos_at(uc, u_curr);
        if (u_ngb_idx < 0) {
            continue;
        }
        const auto np = static_cast<std::size_t>(n_mapped_parent[uc]);
        const auto weight = get_extendable_candidates(uc, u_pos, static_cast<std::size_t>(u_ngb_idx), np, cs);
        if (weight == 0) {
            add_in_failing_set(uc, ancestors[depth].data(), ancestors[depth].size());
            return false;
        }
        if (n_mapped_parent[uc] == 1) {
            queue.insert_to_queue(static_cast<std::int32_t>(uc));
            queue.add_n_inserted(depth);
        }
    }
    return true;
}

void BacktrackContext::update_release_query_ver(std::size_t u_curr, std::size_t depth, const Graph& q) {
    queue.reinsert_to_queue(static_cast<std::int32_t>(u_curr), depth);
    is_mapped[u_curr] = false;
    for (auto i = static_cast<std::size_t>(q.nbr_offset[u_curr]); i < static_cast<std::size_t>(q.nbr_offset[u_curr + 1]); ++i) {
        const auto uc = static_cast<std::size_t>(q.nbr[i].first);
        if (!is_mapped[uc]) {
            --n_mapped_parent[uc];
        }
    }
}

void BacktrackContext::update_failing_set(std::size_t depth) {
    if (is_redundant) {
        ancestors[depth].fill(0);
        return;
    }
    for (std::size_t x = 0; x < failing_set_size; ++x) {
        auto arr = ancestors[depth][x];
        while (arr != 0) {
            const auto bit_pos = static_cast<std::size_t>(std::countr_zero(arr));
            const auto idx = (x << 6) + bit_pos;
            if (idx < kMaxNumVertex) {
                const auto np = static_cast<std::size_t>(n_mapped_parent[idx]);
                for (std::size_t y = 0; y < failing_set_size; ++y) {
                    element[depth].failing_set[y] |= dyanc.get_set_partition(idx, np, y);
                }
            }
            arr &= arr - 1;
        }
        ancestors[depth][x] = 0;
    }
}

bool BacktrackContext::move_up_failing_set(std::size_t depth) {
    if (depth == 0) {
        return false;
    }
    const auto u_id = static_cast<std::size_t>(element[depth].vertex);
    const auto fs_idx = u_id >> 6;
    const auto fs_bit = (1ULL << (u_id & 0x3f));
    if (fs_idx < element[depth + 1].failing_set.size()) {
        if ((element[depth + 1].failing_set[fs_idx] & fs_bit) == 0) {
            for (std::size_t x = 0; x < failing_set_size; ++x) {
                element[depth].failing_set[x] = element[depth + 1].failing_set[x];
            }
            return true;
        }
        for (std::size_t x = 0; x < failing_set_size; ++x) {
            element[depth].failing_set[x] |= element[depth + 1].failing_set[x];
        }
    }
    return false;
}

void BacktrackContext::conflict_class(std::size_t u_curr, std::size_t u_id, std::size_t depth) {
    add_in_failing_set(u_curr, ancestors[depth].data(), ancestors[depth].size());
    add_in_failing_set(u_id, element[depth].failing_set.data(), element[depth].failing_set.size());
    add_in_failing_set(u_id, ancestors[depth].data(), ancestors[depth].size());
}

void BacktrackContext::update_release_candidate(std::size_t depth) {
    const auto v_id = curr_mapping[depth];
    if (v_id >= 0) {
        map_to[static_cast<std::size_t>(v_id)] = -1;
    }
    queue.remove_from_queue(depth);
    queue.clear_n_inserted(depth);
}

void BacktrackContext::backtrack_root_initial(const Graph& q,
                                              const std::vector<CandidateSpace>& cs,
                                              std::size_t root_u) {
    (void)cs;
    const auto nq = static_cast<std::size_t>(q.n_vertex);
    dyanc.clear(root_u);
    queue.clear_queue(nq);
    for (std::size_t i = 0; i < nq; ++i) {
        is_mapped[i] = false;
        n_mapped_parent[i] = 0;
        for (std::size_t j = 0; j < failing_set_size; ++j) {
            ancestors[i][j] = 0;
        }
    }
    is_redundant = false;
    is_mapped[root_u] = true;
    update_ancestors(root_u, q);
    element[0].vertex = static_cast<std::int32_t>(root_u);
    element[0].clear_address();
    element[0].address_pos = -1;
    element[0].address_size = cs[root_u].size;
    for (std::size_t x = 0; x < failing_set_size; ++x) {
        element[0].failing_set[x] = 0ULL;
    }
}

void BacktrackContext::reset_backtrack(std::size_t depth, std::int32_t root_cand, const Graph& q) {
    while (depth > 0) {
        const auto u_id = static_cast<std::size_t>(element[depth].vertex);
        const auto v_id = curr_mapping[depth];
        is_mapped[u_id] = false;
        for (auto i = static_cast<std::size_t>(q.nbr_offset[u_id]); i < static_cast<std::size_t>(q.nbr_offset[u_id + 1]); ++i) {
            const auto uc = static_cast<std::size_t>(q.nbr[i].first);
            --n_mapped_parent[uc];
        }
        queue.clear_n_inserted(depth);
        if (v_id >= 0) {
            map_to[static_cast<std::size_t>(v_id)] = -1;
        }
        element[depth].clear_address();
        --depth;
    }
    if (root_cand >= 0) {
        map_to[static_cast<std::size_t>(root_cand)] = -1;
    }
    queue.clear_queue(static_cast<std::size_t>(q.n_vertex));
    is_redundant = false;
}

bool BacktrackContext::backtrack_once(const Graph& q,
                                      std::vector<CandidateSpace>& cs,
                                      std::size_t root_u,
                                      std::size_t ri) {
    profile::ScopedPhase phase("backtrack_once");
    const auto nq = static_cast<std::size_t>(q.n_vertex);
    backtrack_root_initial(q, cs, root_u);

    const auto root_cand = cs[root_u].candidates[ri];
    map_to[static_cast<std::size_t>(root_cand)] = static_cast<std::int32_t>(root_u);
    curr_mapping[0] = root_cand;
    cand_pos[root_u] = static_cast<std::int32_t>(ri);

    std::size_t depth = 1;
    update_extendable_root(q, cs, ri, root_u);
    bool find_embedding = false;

    while (!find_embedding) {
        if (depth == 0) {
            break;
        }

        if (depth == nq) {
            find_embedding = true;
            newly_confirmed_.clear();
            for (std::size_t internal_depth = 0; internal_depth < nq; ++internal_depth) {
                const auto ui = static_cast<std::size_t>(element[internal_depth].vertex);
                const auto vi_pos = cand_pos[ui];
                if (vi_pos >= 0) {
                    const auto pos = static_cast<std::size_t>(vi_pos);
                    if (pos < static_cast<std::size_t>(cs[ui].size) && cs[ui].marked[pos] == 0) {
                        const auto v = cs[ui].candidates[pos];
                        cs[ui].marked[pos] = 1;
                        ++cs_m[ui];
                        newly_confirmed_.emplace_back(ui, v);
                        sync_mark_confirmed(ui, v, cs);
                    }
                }
            }
            --depth;
            update_release_candidate(depth);
            for (std::size_t x = 0; x < failing_set_size; ++x) {
                element[depth].failing_set[x] = ~0ULL;
            }
            continue;
        }

        if (element[depth].address_size == 0 || element[depth].address_pos < 0) {
            const auto u_curr = static_cast<std::size_t>(queue.pop_from_queue(iec_size, n_mapped_parent, depth));
            element[depth].vertex = static_cast<std::int32_t>(u_curr);
            is_mapped[u_curr] = true;
            update_ancestors(u_curr, q);
            for (std::size_t x = 0; x < failing_set_size; ++x) {
                element[depth].failing_set[x] = 0ULL;
            }
            const auto np = static_cast<std::size_t>(n_mapped_parent[u_curr]);
            if (np == 0 || np > kSmallQueryDegree) {
                element[depth].clear_address();
                is_redundant = false;
                if (depth != 0) {
                    update_failing_set(depth);
                    update_release_query_ver(u_curr, depth, q);
                }
                --depth;
                continue;
            }

            const auto addr_size = get_iec_size(u_curr, np - 1);
            if (addr_size == 0) {
                element[depth].clear_address();
                is_redundant = false;
                if (depth != 0) {
                    update_failing_set(depth);
                    update_release_query_ver(u_curr, depth, q);
                }
                --depth;
                continue;
            }

            std::size_t num_mark0 = 0;
            std::size_t num_mark1 = 0;
            for (std::size_t i = 0; i < static_cast<std::size_t>(addr_size); ++i) {
                const auto pos = get_iec(u_curr, np - 1, i);
                const auto pos_u = static_cast<std::size_t>(pos);
                if (pos_u < static_cast<std::size_t>(cs[u_curr].size)) {
                    if (cs[u_curr].marked[pos_u] == 0) {
                        set_iec(u_curr, np, num_mark0, pos);
                        ++num_mark0;
                    } else if (num_mark1 < mark1_vertex.size()) {
                        mark1_vertex[num_mark1] = pos;
                        ++num_mark1;
                    }
                }
            }

            std::size_t addr_idx = 0;
            for (std::size_t i = 0; i < num_mark0; ++i) {
                element[depth].address[addr_idx++] = get_iec(u_curr, np, i);
            }
            for (std::size_t i = 0; i < num_mark1; ++i) {
                element[depth].address[addr_idx++] = mark1_vertex[i];
            }
            element[depth].address_size = static_cast<std::int32_t>(addr_idx);
            element[depth].address_pos = 0;
        } else {
            const auto u_curr = static_cast<std::size_t>(element[depth].vertex);
            ++element[depth].address_pos;
            if (element[depth].address_pos >= element[depth].address_size || is_redundant) {
                update_failing_set(depth);
                is_redundant = false;
                update_release_query_ver(u_curr, depth, q);
                element[depth].clear_address();
                --depth;
                update_release_candidate(depth);
                is_redundant = move_up_failing_set(depth);
                continue;
            }
        }

        const auto u_curr = static_cast<std::size_t>(element[depth].vertex);
        bool backtrack = false;
        while (true) {
            const auto u_pos_idx = static_cast<std::size_t>(element[depth].address_pos);
            if (u_pos_idx >= static_cast<std::size_t>(element[depth].address_size)) {
                backtrack = true;
                break;
            }
            const auto u_pos = static_cast<std::size_t>(element[depth].address[u_pos_idx]);
            if (u_pos >= static_cast<std::size_t>(cs[u_curr].size)) {
                ++element[depth].address_pos;
                if (element[depth].address_pos >= element[depth].address_size) {
                    backtrack = true;
                    break;
                }
                continue;
            }
            const auto v_id = cs[u_curr].candidates[u_pos];
            if (map_to[static_cast<std::size_t>(v_id)] < 0) {
                curr_mapping[depth] = v_id;
                cand_pos[u_curr] = static_cast<std::int32_t>(u_pos);
                map_to[static_cast<std::size_t>(v_id)] = static_cast<std::int32_t>(u_curr);
                if (!update_extendable_u(q, cs, u_curr, depth, u_pos)) {
                    update_release_candidate(depth);
                    ++element[depth].address_pos;
                    if (element[depth].address_pos >= element[depth].address_size) {
                        backtrack = true;
                        break;
                    }
                    continue;
                }
                break;
            }
            const auto mapped_u = static_cast<std::size_t>(map_to[static_cast<std::size_t>(v_id)]);
            conflict_class(u_curr, mapped_u, depth);
            ++element[depth].address_pos;
            if (element[depth].address_pos >= element[depth].address_size) {
                backtrack = true;
                break;
            }
        }

        if (backtrack) {
            update_failing_set(depth);
            is_redundant = false;
            update_release_query_ver(u_curr, depth, q);
            element[depth].clear_address();
            --depth;
            update_release_candidate(depth);
            is_redundant = move_up_failing_set(depth);
        } else {
            ++depth;
        }
    }

    reset_backtrack(depth, root_cand, q);
    return find_embedding;
}

}
