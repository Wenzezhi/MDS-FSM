#pragma once

#include <cstdint>
#include <vector>

#include "mds/fxhash.hpp"

namespace mds {

struct DfsCode {
    std::int32_t from = 0;
    std::int32_t to = 0;
    std::int32_t from_label = 0;
    std::int32_t edge_label = 0;
    std::int32_t to_label = 0;

    DfsCode() = default;
    DfsCode(std::int32_t from_in,
            std::int32_t to_in,
            std::int32_t from_label_in,
            std::int32_t edge_label_in,
            std::int32_t to_label_in)
        : from(from_in),
          to(to_in),
          from_label(from_label_in),
          edge_label(edge_label_in),
          to_label(to_label_in) {}

    bool operator==(const DfsCode& other) const = default;
    bool operator<(const DfsCode& other) const;
};

std::int32_t build_right_most_path(const std::vector<DfsCode>& dfs_codes,
                                   std::vector<std::int32_t>& right_most_path);

}

namespace mds::fxhash {

inline void append(Hasher& hasher, const mds::DfsCode& value) {
    append(hasher, value.from);
    append(hasher, value.to);
    append(hasher, value.from_label);
    append(hasher, value.edge_label);
    append(hasher, value.to_label);
}

}
