#include "mds/mining/failing_set.hpp"

#include <climits>

namespace mds {

void ExtendableQueue::insert_to_queue(std::int32_t u) {
    extendable_[pos_extendable_] = u;
    ++pos_extendable_;
}

void ExtendableQueue::reinsert_to_queue(std::int32_t u, std::size_t depth) {
    const auto position = static_cast<std::size_t>(positions_[depth]);
    extendable_[pos_extendable_] = extendable_[position];
    ++pos_extendable_;
    extendable_[position] = u;
}

std::int32_t ExtendableQueue::pop_from_queue(const std::vector<std::int32_t>& iec_size,
                                             const std::array<std::int32_t, kMaxQueryVertex>& n_mapped_parent,
                                             std::size_t depth) {
    std::size_t opt_pos = 0;
    std::int64_t opt_weight = LLONG_MAX;

    for (std::size_t i = 0; i < pos_extendable_; ++i) {
        const auto u = static_cast<std::size_t>(extendable_[i]);
        const auto np = n_mapped_parent[u];
        if (np > 0 && static_cast<std::size_t>(np) <= kSmallQueryDegree) {
            const auto idx = u * kSmallQueryDegree + static_cast<std::size_t>(np - 1);
            const auto weight = static_cast<std::int64_t>(iec_size[idx]);
            if (weight < opt_weight) {
                opt_pos = i;
                opt_weight = weight;
            }
        }
    }

    const auto current = extendable_[opt_pos];
    positions_[depth] = static_cast<std::int32_t>(opt_pos);
    extendable_[opt_pos] = extendable_[pos_extendable_ - 1];
    --pos_extendable_;
    return current;
}

void ExtendableQueue::remove_from_queue(std::size_t depth) {
    const auto to_remove = static_cast<std::size_t>(n_inserted_[depth]);
    if (pos_extendable_ >= to_remove) {
        pos_extendable_ -= to_remove;
    }
}

void ExtendableQueue::clear_queue(std::size_t n_query_vertex) {
    for (std::size_t i = 0; i < n_query_vertex; ++i) {
        n_inserted_[i] = 0;
    }
    pos_extendable_ = 0;
}

void ExtendableQueue::clear_n_inserted(std::size_t depth) {
    n_inserted_[depth] = 0;
}

void ExtendableQueue::add_n_inserted(std::size_t depth) {
    ++n_inserted_[depth];
}

DynamicDAGAncestors::DynamicDAGAncestors(std::size_t max_vertex,
                                         std::size_t max_parents,
                                         std::size_t failing_set_size)
    : dag_ancestors_(max_vertex,
                     std::vector<std::vector<std::uint64_t>>(max_parents + 1,
                                                            std::vector<std::uint64_t>(failing_set_size, 0ULL))),
      failing_set_size_(failing_set_size) {}

void DynamicDAGAncestors::clear(std::size_t root) {
    for (std::size_t i = 0; i < failing_set_size_; ++i) {
        dag_ancestors_[root][0][i] = 0ULL;
    }
}

void DynamicDAGAncestors::add_ancestor(std::size_t u,
                                       std::size_t up,
                                       std::size_t u_num_par,
                                       std::size_t up_num_par) {
    if (u_num_par == 1) {
        for (std::size_t i = 0; i < failing_set_size_; ++i) {
            dag_ancestors_[u][1][i] = 0ULL;
        }
    }
    if (u_num_par > 0 && u_num_par < dag_ancestors_[u].size() && up_num_par < dag_ancestors_[up].size()) {
        for (std::size_t x = 0; x < failing_set_size_; ++x) {
            dag_ancestors_[u][u_num_par][x] = dag_ancestors_[u][u_num_par - 1][x];
            dag_ancestors_[u][u_num_par][x] |= dag_ancestors_[up][up_num_par][x];
        }
        dag_ancestors_[u][u_num_par][up >> 6] |= (1ULL << (up & 0x3f));
    }
}

std::uint64_t DynamicDAGAncestors::get_set_partition(std::size_t u,
                                                     std::size_t u_num_par,
                                                     std::size_t y) const {
    if (u < dag_ancestors_.size() && u_num_par < dag_ancestors_[u].size() && y < failing_set_size_) {
        return dag_ancestors_[u][u_num_par][y];
    }
    return 0ULL;
}

void add_in_failing_set(std::size_t u, std::uint64_t* failing_set, std::size_t size) {
    if ((u >> 6) < size) {
        failing_set[u >> 6] |= (1ULL << (u & 0x3f));
    }
}

}
