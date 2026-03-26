#pragma once

#include <cstdint>
#include <optional>
#include <span>
#include <utility>
#include <vector>

#include "mds/candidate_space.hpp"
#include "mds/flat_bitmap.hpp"
#include "mds/graph.hpp"

namespace mds {

class MDSContext {
public:
    MDSContext();

    [[nodiscard]] std::span<const std::size_t> orbit_members(std::size_t oid) const;

    void init(const Graph& q,
              std::vector<CandidateSpace>& cs,
              const std::vector<std::size_t>& orbits,
              std::size_t data_size);

    void mark_confirmed(std::size_t u, std::int32_t v);
    void mark_failed(std::size_t u, std::int32_t v);
    void sync_confirmed_to_cs(std::size_t u, std::int32_t v, std::vector<CandidateSpace>& cs) const;
    void sync_failed_to_cs(std::size_t u, std::int32_t v, std::vector<CandidateSpace>& cs);
    [[nodiscard]] const std::vector<std::pair<std::size_t, std::int32_t>>& get_synced_failed() const;
    void batch_sync_confirmed(const std::vector<CandidateSpace>& cs);
    void update_global_bounds();
    [[nodiscard]] bool is_converged() const;
    [[nodiscard]] bool can_prune(std::int32_t tau) const;
    [[nodiscard]] std::int32_t get_mds() const;
    [[nodiscard]] double success_rate() const;
    [[nodiscard]] std::size_t get_orbit_remaining_count(std::size_t u) const;
    [[nodiscard]] std::optional<std::pair<std::size_t, std::size_t>> select_next(const std::vector<CandidateSpace>& cs);
    void advance_cidx(std::size_t u, std::size_t ci);

    std::vector<std::size_t> orbits;
    std::vector<std::size_t> orbit_size;
    std::size_t num_orbits = 0;
    std::vector<std::int32_t> orbit_label;
    std::vector<std::vector<std::size_t>> label_to_orbits;
    std::size_t num_labels = 0;
    std::vector<std::size_t> active_labels;
    std::int32_t global_lb = 0;
    std::int32_t global_ub = INT32_MAX;
    std::int64_t total_verified = 0;
    std::int64_t total_success = 0;
    double threshold = 0.5;
    std::int32_t tau = 0;

private:
    void update_orbit_intersect_for_failed(std::size_t label_idx, std::size_t oid, std::size_t v_idx);
    [[nodiscard]] std::vector<std::vector<std::size_t>> build_connected_components(
        std::size_t label_idx,
        const std::vector<std::size_t>& orbit_ids) const;
    [[nodiscard]] std::pair<std::int32_t, std::int32_t> compute_component_bounds(
        const std::vector<std::size_t>& component,
        const std::vector<std::size_t>& orbit_ids,
        std::size_t label_idx);
    [[nodiscard]] std::pair<std::int32_t, std::int32_t> compute_class_bounds(std::size_t label_idx);
    [[nodiscard]] std::optional<std::pair<std::size_t, std::size_t>> select_for_raise_lb(
        const std::vector<CandidateSpace>& cs);
    [[nodiscard]] std::optional<std::pair<std::size_t, std::size_t>> select_for_lower_ub(
        const std::vector<CandidateSpace>& cs);

    std::vector<std::size_t> orbit_members_data_;
    std::vector<std::size_t> orbit_members_offset_;
    std::size_t bitmap_words_ = 0;
    FlatBitmap orbit_confirmed_;
    FlatBitmap orbit_remaining_;
    std::vector<std::size_t> orbit_confirmed_count_;
    std::vector<std::size_t> orbit_remaining_count_;
    std::vector<std::int32_t> orbit_lb_;
    std::vector<std::int32_t> orbit_ub_;
    std::vector<std::int32_t> class_lb_;
    std::vector<std::int32_t> class_ub_;
    std::vector<std::uint8_t> class_dirty_;
    std::vector<std::uint64_t> temp_union_;
    std::vector<std::uint64_t> temp_union2_;
    std::vector<std::size_t> cs_cidx_;
    std::vector<std::pair<std::size_t, std::int32_t>> sync_buffer_;
    std::vector<std::vector<std::size_t>> class_subset_denom_;
    std::vector<std::vector<std::uint8_t>> orbit_intersect_;
    std::vector<std::size_t> orbit_index_in_label_;
    std::vector<std::uint64_t> temp_intersect_;
    std::vector<std::uint8_t> orbit_has_unverified_;
    std::vector<std::size_t> orbit_unverified_count_;
    std::vector<std::size_t> vertex_unverified_count_;
};

}
