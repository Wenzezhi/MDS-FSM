#include "mds/candidate_space.hpp"

namespace mds {

void CandidateSpace::init(std::size_t max_candidates, std::size_t query_degree, std::size_t data_n_vertex) {
    size = 0;
    candidates.assign(max_candidates, 0);
    marked.assign(max_candidates, 0);
    reset_vertex_to_pos(data_n_vertex);
    adjacent = FlatMatrix<std::vector<std::int32_t>>(query_degree, max_candidates);
    nbr_cnt = FlatMatrix<std::int32_t>(max_candidates, query_degree);
    nbr_set_cnt = FlatMatrix<std::int8_t>(0, 0);
    nbr_safety = FlatMatrix<std::int32_t>(max_candidates, query_degree);
}

void CandidateSpace::init_nbr_set_cnt(std::size_t max_candidates, std::size_t max_data_degree) {
    nbr_set_cnt = FlatMatrix<std::int8_t>(max_candidates, max_data_degree);
}

void CandidateSpace::reset_vertex_to_pos(std::size_t data_n_vertex) {
    if (vertex_to_pos.size() != data_n_vertex) {
        vertex_to_pos.assign(data_n_vertex, 0);
        vertex_to_pos_stamp.assign(data_n_vertex, 0);
        vertex_to_pos_epoch = 1;
        return;
    }

    ++vertex_to_pos_epoch;
    if (vertex_to_pos_epoch == 0) {
        std::fill(vertex_to_pos_stamp.begin(), vertex_to_pos_stamp.end(), 0U);
        vertex_to_pos_epoch = 1;
    }
}

void CandidateSpace::set_vertex_pos(std::size_t vertex, std::int32_t pos) {
    vertex_to_pos[vertex] = pos;
    vertex_to_pos_stamp[vertex] = vertex_to_pos_epoch;
}

std::int32_t CandidateSpace::get_vertex_pos(std::size_t vertex) const {
    return vertex < vertex_to_pos.size() && vertex_to_pos_stamp[vertex] == vertex_to_pos_epoch ? vertex_to_pos[vertex]
                                                                                                : -1;
}

void CandidateSpace::clear_vertex_pos(std::size_t vertex) {
    if (vertex < vertex_to_pos_stamp.size()) {
        vertex_to_pos_stamp[vertex] = 0;
    }
}

}
