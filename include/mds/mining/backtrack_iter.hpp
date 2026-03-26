#pragma once

#include <array>
#include <cstdint>
#include <span>
#include <utility>
#include <vector>

#include "mds/candidate_space.hpp"
#include "mds/graph.hpp"
#include "mds/mining/failing_set.hpp"
#include "mds/types.hpp"

namespace mds {

class BacktrackContext {
public:
    BacktrackContext(std::size_t max_data_vertex, std::size_t max_candidate, std::size_t max_degree);

    void set_orbits(const std::vector<std::size_t>& orbits);
    void sync_mark_confirmed(std::size_t u, std::int32_t v, std::vector<CandidateSpace>& cs);
    [[nodiscard]] std::vector<std::pair<std::size_t, std::int32_t>> sync_mark_failed(
        std::size_t u,
        std::int32_t v,
        std::vector<CandidateSpace>& cs);
    [[nodiscard]] const std::vector<std::pair<std::size_t, std::int32_t>>& get_newly_confirmed() const;
    void clear_newly_confirmed();
    [[nodiscard]] bool backtrack_once(const Graph& q,
                                      std::vector<CandidateSpace>& cs,
                                      std::size_t root_u,
                                      std::size_t ri);

    std::vector<StackFrame> element;
    ExtendableQueue queue;
    std::vector<std::array<std::uint64_t, 1>> ancestors;
    DynamicDAGAncestors dyanc;
    std::size_t failing_set_size = 0;
    bool is_redundant = false;
    std::vector<std::int32_t> mark1_vertex;
    std::vector<std::int32_t> map_to;
    std::array<bool, kMaxQueryVertex> is_mapped{};
    std::array<std::int32_t, kMaxQueryVertex> curr_mapping{};
    std::array<std::int32_t, kMaxQueryVertex> cand_pos{};
    std::array<std::int32_t, kMaxQueryVertex> n_mapped_parent{};
    std::vector<std::int32_t> iec_data;
    std::vector<std::int32_t> iec_size;
    std::size_t iec_capacity = 0;
    std::array<std::int32_t, kMaxQueryVertex> cs_m{};
    std::array<std::int32_t, kMaxQueryVertex> cs_c{};
    std::array<std::int32_t, kMaxQueryVertex> cs_cidx{};
    std::vector<std::size_t> orbits;

private:
    [[nodiscard]] std::size_t iec_idx(std::size_t u, std::size_t p, std::size_t i) const;
    [[nodiscard]] static std::size_t iec_size_idx(std::size_t u, std::size_t p);
    [[nodiscard]] std::int32_t get_iec_size(std::size_t u, std::size_t p) const;
    void set_iec_size(std::size_t u, std::size_t p, std::int32_t value);
    [[nodiscard]] std::int32_t get_iec(std::size_t u, std::size_t p, std::size_t i) const;
    void set_iec(std::size_t u, std::size_t p, std::size_t i, std::int32_t value);
    [[nodiscard]] std::int64_t get_extendable_candidates(std::size_t u,
                                                         std::size_t cand_parent_idx,
                                                         std::size_t u_ngb_idx,
                                                         std::size_t num_parents_mapped,
                                                         const std::vector<CandidateSpace>& cs);
    void update_ancestors(std::size_t u_curr, const Graph& q);
    void update_extendable_root(const Graph& q,
                                const std::vector<CandidateSpace>& cs,
                                std::size_t ri,
                                std::size_t root_u);
    [[nodiscard]] bool update_extendable_u(const Graph& q,
                                           const std::vector<CandidateSpace>& cs,
                                           std::size_t u_curr,
                                           std::size_t depth,
                                           std::size_t u_pos);
    void update_release_query_ver(std::size_t u_curr, std::size_t depth, const Graph& q);
    void update_failing_set(std::size_t depth);
    [[nodiscard]] bool move_up_failing_set(std::size_t depth);
    void conflict_class(std::size_t u_curr, std::size_t u_id, std::size_t depth);
    void update_release_candidate(std::size_t depth);
    void backtrack_root_initial(const Graph& q, const std::vector<CandidateSpace>& cs, std::size_t root_u);
    void reset_backtrack(std::size_t depth, std::int32_t root_cand, const Graph& q);
    [[nodiscard]] std::span<const std::size_t> orbit_members_span(std::size_t oid) const;

    std::vector<std::pair<std::size_t, std::int32_t>> newly_confirmed_;
    std::vector<std::size_t> orbit_members_data_;
    std::vector<std::size_t> orbit_members_offset_;
};

}
