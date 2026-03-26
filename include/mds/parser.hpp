#pragma once

#include <cstdint>
#include <string>
#include <tuple>
#include <unordered_map>

#include "mds/fxhash.hpp"
#include "mds/graph.hpp"

namespace mds {

using LabelMap = std::unordered_map<std::string, std::int32_t, fxhash::Hash<std::string>>;
using LabelInverseMap = std::unordered_map<std::int32_t, std::string, fxhash::Hash<std::int32_t>>;
using EdgeLabelMap = std::unordered_map<std::int32_t, std::int32_t, fxhash::Hash<std::int32_t>>;
using EdgeListMap = std::unordered_map<std::tuple<std::int32_t, std::int32_t, std::int32_t>,
                                       std::int32_t,
                                       fxhash::Hash<std::tuple<std::int32_t, std::int32_t, std::int32_t>>>;

struct ParseResult {
    Graph graph;
    LabelMap label_map;
    LabelInverseMap label_map_inverse;
    EdgeLabelMap edge_label_map;
    EdgeLabelMap edge_label_map_inverse;
    EdgeListMap edge_list;
};

ParseResult read_gfu_format(const std::string& file_path);

}
