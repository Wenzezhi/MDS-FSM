#pragma once

#include <cstdint>
#include <vector>

#include "mds/candidate_space.hpp"
#include "mds/graph.hpp"

namespace mds {

class MiningContext;

std::vector<std::size_t> compute_orbits(const Graph& q);
std::int32_t compute_mds_ub(const Graph& q,
                           const std::vector<CandidateSpace>& cs,
                           const std::vector<std::size_t>& orbits,
                           MiningContext& ctx);
std::int32_t compute_mds_from_cs(const Graph& q,
                                 const std::vector<CandidateSpace>& cs,
                                 const std::vector<std::size_t>& orbits,
                                 MiningContext& ctx);

}
