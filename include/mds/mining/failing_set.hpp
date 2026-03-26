#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "mds/types.hpp"

namespace mds {

inline constexpr std::size_t kMaxStackCandidates = 65536;
inline constexpr std::size_t kMaxFailingSetSize = 8;

struct StackFrame {
    std::array<std::int32_t, kMaxStackCandidates> address{};
    std::int32_t address_size = 0;
    std::int32_t address_pos = -1;
    std::int32_t vertex = -1;
    std::array<std::uint64_t, kMaxFailingSetSize> failing_set{};

    explicit StackFrame(std::size_t failing_set_size) {
        (void)failing_set_size;
    }

    void clear_address() {
        address_size = 0;
        address_pos = -1;
    }
};

class ExtendableQueue {
public:
    ExtendableQueue() = default;

    void insert_to_queue(std::int32_t u);
    void reinsert_to_queue(std::int32_t u, std::size_t depth);
    std::int32_t pop_from_queue(const std::vector<std::int32_t>& iec_size,
                                const std::array<std::int32_t, kMaxQueryVertex>& n_mapped_parent,
                                std::size_t depth);
    void remove_from_queue(std::size_t depth);
    void clear_queue(std::size_t n_query_vertex);
    void clear_n_inserted(std::size_t depth);
    void add_n_inserted(std::size_t depth);

private:
    std::size_t pos_extendable_ = 0;
    std::array<std::int32_t, kMaxNumVertex> extendable_{};
    std::array<std::int32_t, kMaxNumVertex> positions_{};
    std::array<std::int32_t, kMaxNumVertex> n_inserted_{};
};

class DynamicDAGAncestors {
public:
    DynamicDAGAncestors(std::size_t max_vertex, std::size_t max_parents, std::size_t failing_set_size);

    void clear(std::size_t root);
    void add_ancestor(std::size_t u, std::size_t up, std::size_t u_num_par, std::size_t up_num_par);
    [[nodiscard]] std::uint64_t get_set_partition(std::size_t u, std::size_t u_num_par, std::size_t y) const;

private:
    std::vector<std::vector<std::vector<std::uint64_t>>> dag_ancestors_;
    std::size_t failing_set_size_ = 0;
};

void add_in_failing_set(std::size_t u, std::uint64_t* failing_set, std::size_t size);

}
