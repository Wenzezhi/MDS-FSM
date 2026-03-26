#pragma once

#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "mds/candidate_space.hpp"
#include "mds/edge_image_bitmaps.hpp"
#include "mds/graph.hpp"
#include "mds/minimality_cache.hpp"
#include "mds/object_pool.hpp"
#include "mds/parser.hpp"
#include "mds/mining/backtrack_iter.hpp"
#include "mds/mining/heap.hpp"

namespace mds {

using FastIntSet = std::unordered_set<std::int32_t, fxhash::Hash<std::int32_t>>;
using OrbitCache = std::unordered_map<std::uint64_t, std::vector<std::size_t>, fxhash::Hash<std::uint64_t>>;

class MiningContext {
public:
    MiningContext();

    void load_graph(const std::string& path);
    void process_data_graph();
    void compute_frequent_edges(std::size_t k);
    void frequent_mining(std::size_t k, std::size_t tth);
    void print_answer() const;

    [[nodiscard]] Graph create_subgraph_edge(std::int32_t l1, std::int32_t el, std::int32_t l2) const;
    [[nodiscard]] Graph create_subgraph(const Graph& p,
                                        bool fwd,
                                        std::int32_t new_lbl,
                                        std::int32_t el,
                                        std::int32_t fst,
                                        std::int32_t snd) const;
    [[nodiscard]] bool is_in_freq_edge(std::int32_t l1, std::int32_t el, std::int32_t l2) const;
    [[nodiscard]] bool is_min(const Graph& g);
    [[nodiscard]] std::vector<CandidateSpace> alloc_cs(const Graph& q) const;
    [[nodiscard]] bool filter_count(const Graph& q) const;
    [[nodiscard]] bool build_cs_edge(const Graph& q, std::vector<CandidateSpace>& cs);
    [[nodiscard]] std::int32_t compute_ub(const Graph& q, const std::vector<CandidateSpace>& cs) const;
    [[nodiscard]] std::optional<std::pair<std::vector<CandidateSpace>, std::int32_t>> build_cs_from_edge(Graph& q);
    [[nodiscard]] std::optional<std::vector<CandidateSpace>> copy_cs_fwd(const Graph& q,
                                                                         const Graph& p,
                                                                         const std::vector<CandidateSpace>& pcs);
    [[nodiscard]] std::optional<std::vector<CandidateSpace>> copy_cs_bwd(const Graph& q,
                                                                         const Graph& p,
                                                                         const std::vector<CandidateSpace>& pcs);
    [[nodiscard]] std::optional<std::pair<std::vector<CandidateSpace>, std::int32_t>> build_cs_fwd(
        Graph& q,
        const Graph& p,
        const std::vector<CandidateSpace>& pcs);
    [[nodiscard]] std::optional<std::pair<std::vector<CandidateSpace>, std::int32_t>> build_cs_bwd(
        Graph& q,
        const Graph& p,
        const std::vector<CandidateSpace>& pcs);

    void build_nbr_cnt(const Graph& q, std::vector<CandidateSpace>& cs) const;
    [[nodiscard]] bool cs_node_filtering(const Graph& q, std::vector<CandidateSpace>& cs) const;
    [[nodiscard]] bool cs_node_filtering2(const Graph& q,
                                          std::vector<CandidateSpace>& cs,
                                          std::vector<std::int32_t>& cs_c,
                                          std::size_t iu,
                                          std::int32_t iv,
                                          std::vector<std::pair<std::size_t, std::int32_t>>& newly_failed) const;

    [[nodiscard]] std::string format_pattern(const Graph& g) const;
    [[nodiscard]] const std::vector<std::size_t>& get_or_compute_orbits(const Graph& q);
    [[nodiscard]] std::uint64_t extend_edge(CandidateHeap& gs, const Graph& p, std::int32_t parent_mds);
    [[nodiscard]] std::uint64_t extend(CandidateHeap& gs,
                                       const Graph& p,
                                       const std::vector<CandidateSpace>& pcs,
                                       std::int32_t mni);
    [[nodiscard]] bool prepare_mni(const Graph& q, std::vector<CandidateSpace>& cs);
    [[nodiscard]] std::int32_t compute_mni(const Graph& q, std::vector<CandidateSpace>& cs);

    Graph data_graph;
    LabelMap label_map;
    LabelInverseMap label_map_inverse;
    EdgeLabelMap edge_label_map;
    EdgeLabelMap edge_label_map_inverse;
    EdgeListMap edge_list;
    std::size_t n_unique_label = 0;
    std::size_t n_unique_edge_label = 0;
    std::size_t nlf_size = 0;
    std::size_t max_num_candidate = 0;
    std::size_t max_degree = 0;
    std::size_t max_num_data_vertex = 0;
    std::int32_t tau = 0;
    bool use_filter = false;
    std::vector<std::tuple<std::int32_t, std::int32_t, std::int32_t, std::int32_t>> freq_edge_set;
    std::vector<std::tuple<std::int32_t, std::int32_t, std::int32_t>> freq_edge;
    std::vector<std::int32_t> freq_edge_offset;
    std::vector<std::pair<FastIntSet, FastIntSet>> edge_img_sets;
    std::optional<EdgeImageBitmaps> edge_img_bitmaps;
    bool use_bitmap_optimization = true;
    std::vector<std::int32_t> inv_nbr;
    std::vector<std::uint8_t> is_visited;
    std::vector<std::int32_t> array_to_clean;
    std::size_t to_clean_index = 0;
    std::uint64_t patterns_verified = 0;
    std::uint64_t candidates_generated = 0;
    std::uint64_t candidates_added_to_queue = 0;
    std::vector<AnswerElement> answer_set;
    std::unique_ptr<BacktrackContext> bt_ctx;
    OrbitCache orbit_cache;
    MinimalityCache minimality_cache;
    ObjectPool<std::vector<std::int32_t>> temp_vec_pool;
    ObjectPool<FastIntSet> temp_set_pool;
    std::vector<std::int32_t> cs_c_buffer;
    mutable std::vector<std::pair<std::size_t, std::int32_t>> filter_deadnodes_buffer;
    mutable std::vector<std::int32_t> filter_cs_bc_buffer;

private:
    void compute_core_numbers();
    void build_inv_nbr();
    void allocate_structures();
    void compute_frequent_edges_bitmap(std::size_t k);
    void compute_frequent_edges_hashset(std::size_t k);
    void build_freq_edge_structures(std::size_t k);
    [[nodiscard]] bool is_min_uncached(const Graph& g);
    [[nodiscard]] bool build_cmp_dfs(const Graph& g,
                                     std::vector<std::int32_t>& rn,
                                     std::vector<std::int32_t>& ri,
                                     std::int32_t cn,
                                     std::size_t dl,
                                     std::vector<std::int32_t>& stk) const;
    void remove_node(std::size_t u, std::int32_t v, std::vector<CandidateSpace>& cs, const Graph& q) const;
};

}
