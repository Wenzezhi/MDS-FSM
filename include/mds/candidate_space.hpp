#pragma once

#include <cstdint>
#include <vector>

#include "mds/flat_matrix.hpp"

namespace mds {

struct CandidateSpace {
    std::int32_t size = 0;
    std::vector<std::int32_t> candidates;
    std::vector<std::int8_t> marked;
    std::vector<std::int32_t> vertex_to_pos;
    std::vector<std::uint32_t> vertex_to_pos_stamp;
    std::uint32_t vertex_to_pos_epoch = 1;
    FlatMatrix<std::vector<std::int32_t>> adjacent;
    FlatMatrix<std::int32_t> nbr_cnt;
    FlatMatrix<std::int8_t> nbr_set_cnt;
    FlatMatrix<std::int32_t> nbr_safety;

    void init(std::size_t max_candidates, std::size_t query_degree, std::size_t data_n_vertex);
    void init_nbr_set_cnt(std::size_t max_candidates, std::size_t max_data_degree);
    void reset_vertex_to_pos(std::size_t data_n_vertex);
    void set_vertex_pos(std::size_t vertex, std::int32_t pos);
    [[nodiscard]] std::int32_t get_vertex_pos(std::size_t vertex) const;
    void clear_vertex_pos(std::size_t vertex);

    [[nodiscard]] std::size_t nbr_cnt_len() const { return nbr_cnt.num_rows(); }
    [[nodiscard]] std::size_t nbr_cnt_row_len() const { return nbr_cnt.row_capacity(); }
    [[nodiscard]] std::size_t nbr_set_cnt_len() const { return nbr_set_cnt.num_rows(); }
    [[nodiscard]] std::size_t nbr_set_cnt_row_len() const { return nbr_set_cnt.row_capacity(); }
    [[nodiscard]] std::size_t nbr_safety_len() const { return nbr_safety.num_rows(); }
    [[nodiscard]] std::size_t nbr_safety_row_len() const { return nbr_safety.row_capacity(); }
};

}
