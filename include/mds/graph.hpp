#pragma once

#include <cstdint>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include "mds/dfs_code.hpp"
#include "mds/fxhash.hpp"
#include "mds/types.hpp"

namespace mds {

using LabelFrequencyMap = std::unordered_map<std::int32_t, std::int32_t, fxhash::Hash<std::int32_t>>;

struct Graph {
    NumberType n_vertex = 0;
    NumberType n_edge = 0;
    std::int32_t n_label = 0;
    std::int32_t max_degree = 0;
    std::int32_t max_edge_label = 0;

    std::vector<std::int32_t> label;
    std::vector<std::pair<std::int32_t, std::int32_t>> nbr;
    std::vector<NumberType> nbr_offset;
    std::vector<std::int32_t> degree;
    std::vector<std::int32_t> vertex;
    std::vector<std::int32_t> vertex_offset;
    std::vector<std::int32_t> max_nbr_degree;
    std::vector<std::int32_t> core;

    std::vector<std::uint64_t> nlf;
    LabelFrequencyMap label_frequency;

    std::vector<std::int32_t> right_most_path;
    std::int32_t rmp_size = 0;

    std::vector<DfsCode> dfs_code;

    std::vector<std::int32_t> nbr_to_pos_data;
    std::vector<std::int32_t> comp_label_idx_data;
    std::vector<std::int32_t> nbr_label_data;
    std::vector<std::int32_t> nbr_comp_idx_data;
    std::vector<std::int32_t> nbr_reverse_pos_data;
    std::vector<std::int32_t> nbr_reverse_comp_idx_data;
    std::vector<std::pair<NumberType, NumberType>> label_nbr_offset_table;
    std::size_t label_nbr_offset_stride = 0;

    void init_vertex(NumberType n);
    void init_nbr_to_pos();
    void init_edge(NumberType n);
    void sort_neighbors(std::size_t n_unique_label);
    void compute_nlf(std::size_t nlf_size, std::size_t n_unique_label);
    void compute_nbr_safety();
    [[nodiscard]] std::pair<NumberType, NumberType> lookup_label_nbr_offset(std::int32_t vertex_id,
                                                                             std::int32_t label_id) const;
    [[nodiscard]] std::int32_t nbr_to_pos_at(std::size_t u, std::size_t v) const;
    void set_nbr_to_pos(std::size_t u, std::size_t v, std::int32_t value);
    [[nodiscard]] std::int32_t comp_label_idx_at(std::size_t u, std::size_t v) const;
    void set_comp_label_idx(std::size_t u, std::size_t v, std::int32_t value);

private:
    [[nodiscard]] std::size_t pair_matrix_index(std::size_t u, std::size_t v) const;
};

std::uint64_t hash_pattern(const std::vector<DfsCode>& dfs_code);

}
