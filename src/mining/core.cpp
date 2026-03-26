#include "mds/mining/core.hpp"

#include <algorithm>
#include <cstdio>
#include <iostream>

#include "mds/mining/delayed_queue.hpp"
#include "mds/mining/heap.hpp"
#include "mds/profiler.hpp"

namespace mds {

MiningContext::MiningContext()
    : minimality_cache(MinimalityCache::with_default_size()),
      temp_vec_pool(32),
      temp_set_pool(16) {}

void MiningContext::load_graph(const std::string& path) {
    auto parsed = read_gfu_format(path);
    data_graph = std::move(parsed.graph);
    label_map = std::move(parsed.label_map);
    label_map_inverse = std::move(parsed.label_map_inverse);
    edge_label_map = std::move(parsed.edge_label_map);
    edge_label_map_inverse = std::move(parsed.edge_label_map_inverse);
    edge_list = std::move(parsed.edge_list);
}

void MiningContext::process_data_graph() {
    n_unique_label = label_map.size();
    n_unique_edge_label = edge_label_map.size();
    use_filter = static_cast<double>(data_graph.n_edge) / static_cast<double>(data_graph.n_vertex) >= 1.5;
    nlf_size = n_unique_label == 0 ? 1 : (kBitsPerLabel * n_unique_label - 1) / 64 + 1;
    max_degree = static_cast<std::size_t>(data_graph.max_degree);
    max_num_data_vertex = static_cast<std::size_t>(data_graph.n_vertex);

    std::int32_t max_val = 0;
    for (const auto& [label, count] : data_graph.label_frequency) {
        (void)label;
        max_val = std::max(max_val, count);
    }
    max_num_candidate = static_cast<std::size_t>(max_val);

    compute_core_numbers();
    data_graph.sort_neighbors(n_unique_label);
    data_graph.compute_nlf(nlf_size, n_unique_label);
    build_inv_nbr();
    allocate_structures();
}

void MiningContext::compute_core_numbers() {
    const auto n = static_cast<std::size_t>(data_graph.n_vertex);
    const auto max_d = static_cast<std::size_t>(data_graph.max_degree);
    std::vector<std::int32_t> bin(max_d + 1, 0);
    std::vector<std::size_t> pos(n, 0);
    std::vector<std::size_t> vert(n, 0);

    for (std::size_t i = 0; i < n; ++i) {
        ++bin[static_cast<std::size_t>(data_graph.core[i])];
    }
    std::int32_t start = 0;
    for (std::size_t d = 0; d <= max_d; ++d) {
        const auto num = bin[d];
        bin[d] = start;
        start += num;
    }
    for (std::size_t i = 0; i < n; ++i) {
        pos[i] = static_cast<std::size_t>(bin[static_cast<std::size_t>(data_graph.core[i])]);
        vert[pos[i]] = i;
        ++bin[static_cast<std::size_t>(data_graph.core[i])];
    }
    for (std::size_t d = max_d; d > 0; --d) {
        bin[d] = bin[d - 1];
    }
    bin[0] = 0;

    for (std::size_t i = 0; i < n; ++i) {
        const auto v = vert[i];
        for (auto j = static_cast<std::size_t>(data_graph.nbr_offset[v]);
             j < static_cast<std::size_t>(data_graph.nbr_offset[v + 1]);
             ++j) {
            const auto u = static_cast<std::size_t>(data_graph.nbr[j].first);
            if (data_graph.core[u] > data_graph.core[v]) {
                const auto du = static_cast<std::size_t>(data_graph.core[u]);
                const auto pu = pos[u];
                const auto pw = static_cast<std::size_t>(bin[du]);
                const auto w = vert[static_cast<std::size_t>(bin[du])];
                if (u != w) {
                    pos[u] = pw;
                    pos[w] = pu;
                    vert[pu] = w;
                    vert[pw] = u;
                }
                ++bin[du];
                --data_graph.core[u];
            }
        }
    }
}

void MiningContext::build_inv_nbr() {
    const auto n = max_num_data_vertex;
    inv_nbr.assign(data_graph.nbr.size(), 0);
    std::vector<std::int32_t> curr_idx(n, 0);
    for (std::size_t i = 0; i < n; ++i) {
        const auto v = static_cast<std::size_t>(data_graph.vertex[i]);
        for (auto j = static_cast<std::size_t>(data_graph.nbr_offset[v]);
             j < static_cast<std::size_t>(data_graph.nbr_offset[v + 1]);
             ++j) {
            const auto vn = static_cast<std::size_t>(data_graph.nbr[j].first);
            inv_nbr[j] = curr_idx[vn];
            ++curr_idx[vn];
        }
    }
}

void MiningContext::allocate_structures() {
    const auto n = max_num_data_vertex;
    is_visited.assign(n, 0);
    array_to_clean.assign(n, 0);
    filter_deadnodes_buffer.clear();
    filter_deadnodes_buffer.reserve(std::max<std::size_t>(max_num_candidate, 64));
    filter_cs_bc_buffer.assign(64, 0);
    bt_ctx = std::make_unique<BacktrackContext>(n, max_num_candidate, max_degree);
}

void MiningContext::compute_frequent_edges(std::size_t k) {
    profile::ScopedPhase phase("frequent_edges");
    if (use_bitmap_optimization) {
        compute_frequent_edges_bitmap(k);
    } else {
        compute_frequent_edges_hashset(k);
    }
}

void MiningContext::compute_frequent_edges_bitmap(std::size_t k) {
    const auto edge_list_num = edge_list.size();
    const auto n_vertex = static_cast<std::size_t>(data_graph.n_vertex);
    EdgeImageBitmaps bitmaps(edge_list_num, n_vertex);

    for (std::size_t i = 0; i < n_vertex; ++i) {
        for (auto j = static_cast<std::size_t>(data_graph.nbr_offset[i]);
             j < static_cast<std::size_t>(data_graph.nbr_offset[i + 1]);
             ++j) {
            const auto v = static_cast<std::size_t>(data_graph.nbr[j].first);
            const auto el = data_graph.nbr[j].second;
            const auto vl = data_graph.label[v];
            if (data_graph.label[i] <= vl) {
                const auto key = std::make_tuple(data_graph.label[i], vl, el);
                const auto it = edge_list.find(key);
                if (it != edge_list.end()) {
                    const auto eid = static_cast<std::size_t>(it->second);
                    bitmaps.insert(eid, 0, i);
                    bitmaps.insert(eid, 1, v);
                }
            }
        }
    }

    for (const auto& [key, eid] : edge_list) {
        const auto mds = bitmaps.compute_mds(static_cast<std::size_t>(eid));
        if (mds > 0) {
            freq_edge_set.emplace_back(mds, std::get<0>(key), std::get<2>(key), std::get<1>(key));
        }
    }

    edge_img_bitmaps = std::move(bitmaps);
    std::stable_sort(freq_edge_set.begin(), freq_edge_set.end(), [](const auto& a, const auto& b) {
        return std::get<0>(b) < std::get<0>(a);
    });
    build_freq_edge_structures(k);
}

void MiningContext::compute_frequent_edges_hashset(std::size_t k) {
    edge_img_sets.assign(edge_list.size(), {});
    for (std::size_t i = 0; i < static_cast<std::size_t>(data_graph.n_vertex); ++i) {
        for (auto j = static_cast<std::size_t>(data_graph.nbr_offset[i]);
             j < static_cast<std::size_t>(data_graph.nbr_offset[i + 1]);
             ++j) {
            const auto v = static_cast<std::size_t>(data_graph.nbr[j].first);
            const auto el = data_graph.nbr[j].second;
            const auto vl = data_graph.label[v];
            if (data_graph.label[i] <= vl) {
                const auto it = edge_list.find(std::make_tuple(data_graph.label[i], vl, el));
                if (it != edge_list.end()) {
                    edge_img_sets[static_cast<std::size_t>(it->second)].first.insert(static_cast<std::int32_t>(i));
                    edge_img_sets[static_cast<std::size_t>(it->second)].second.insert(static_cast<std::int32_t>(v));
                }
            }
        }
    }

    for (const auto& [key, eid] : edge_list) {
        const auto& img1 = edge_img_sets[static_cast<std::size_t>(eid)].first;
        const auto& img2 = edge_img_sets[static_cast<std::size_t>(eid)].second;
        FastIntSet union_set = img1;
        union_set.insert(img2.begin(), img2.end());
        const auto mds = static_cast<std::int32_t>(
            std::min({union_set.size() / 2, img1.size(), img2.size()}));
        if (mds > 0) {
            freq_edge_set.emplace_back(mds, std::get<0>(key), std::get<2>(key), std::get<1>(key));
        }
    }

    std::stable_sort(freq_edge_set.begin(), freq_edge_set.end(), [](const auto& a, const auto& b) {
        return std::get<0>(b) < std::get<0>(a);
    });
    build_freq_edge_structures(k);
}

void MiningContext::build_freq_edge_structures(std::size_t k) {
    if (freq_edge_set.size() > k) {
        const auto k_mds = std::get<0>(freq_edge_set[k - 1]);
        auto end_pos = k;
        while (end_pos < freq_edge_set.size() && std::get<0>(freq_edge_set[end_pos]) == k_mds) {
            ++end_pos;
        }
        freq_edge_set.resize(end_pos);
    }

    std::vector<std::int32_t> label_freq(n_unique_label, 0);
    for (const auto& [m, x, l, y] : freq_edge_set) {
        (void)m;
        (void)l;
        ++label_freq[static_cast<std::size_t>(x)];
        if (x != y) {
            ++label_freq[static_cast<std::size_t>(y)];
        }
    }

    freq_edge_offset.assign(n_unique_label + 1, 0);
    std::int32_t temp = 0;
    for (std::size_t i = 0; i < n_unique_label; ++i) {
        freq_edge_offset[i] = temp;
        temp += label_freq[i];
    }
    freq_edge_offset[n_unique_label] = temp;

    freq_edge.assign(static_cast<std::size_t>(temp), {0, 0, 0});
    for (const auto& [m, x, l, y] : freq_edge_set) {
        auto pos = static_cast<std::size_t>(freq_edge_offset[static_cast<std::size_t>(x)]);
        freq_edge[pos] = {m, l, y};
        ++freq_edge_offset[static_cast<std::size_t>(x)];
        if (x != y) {
            pos = static_cast<std::size_t>(freq_edge_offset[static_cast<std::size_t>(y)]);
            freq_edge[pos] = {m, l, x};
            ++freq_edge_offset[static_cast<std::size_t>(y)];
        }
    }
    for (std::size_t i = n_unique_label; i > 0; --i) {
        freq_edge_offset[i] = freq_edge_offset[i - 1];
    }
    freq_edge_offset[0] = 0;
}

void MiningContext::frequent_mining(std::size_t k, std::size_t tth) {
    candidates_generated = 0;
    candidates_added_to_queue = 0;
    patterns_verified = 0;

    compute_frequent_edges(k);

    CandidateHeap cand_graph_set;
    AnswerHeap answer_heap;
    std::vector<std::pair<Graph, std::int32_t>> freq_edge_graphs;
    freq_edge_graphs.reserve(k);
    DelayedExtensionQueue delayed_queue(tth);

    for (std::size_t i = 0; i < std::min(k, freq_edge_set.size()); ++i) {
        const auto [mds, l1, el, l2] = freq_edge_set[i];
        auto g = create_subgraph_edge(l1, el, l2);
        answer_heap.push(AnswerElement{mds, g});
        freq_edge_graphs.emplace_back(std::move(g), mds);
    }

    tau = answer_heap.size() == k ? answer_heap.top().mds : 0;

    for (const auto& [edge_graph, mds] : freq_edge_graphs) {
        const auto tau_tth = get_tth_largest_ub(cand_graph_set, tth);
        if (mds > tau_tth) {
            profile::ScopedPhase phase("extend_total");
            candidates_added_to_queue += extend_edge(cand_graph_set, edge_graph, mds);
            delayed_queue.record_immediate();
        } else {
            delayed_queue.push(mds, edge_graph);
        }
    }

    while (!cand_graph_set.empty()) {
        auto helt = pop_heap_element(cand_graph_set);
        if (helt.ub <= tau) {
            break;
        }

        auto payload = std::move(helt.payload);
        auto graph = std::move(payload->graph);
        auto cand_space = std::move(payload->cand_space);
        if (!prepare_mni(graph, cand_space)) {
            continue;
        }

        ++patterns_verified;
        const auto new_mds = compute_mni(graph, cand_space);
        if (new_mds > tau) {
            if (answer_heap.size() < k) {
                answer_heap.push(AnswerElement{new_mds, graph});
            } else {
                answer_heap.pop();
                answer_heap.push(AnswerElement{new_mds, graph});
            }
            if (answer_heap.size() == k) {
                tau = answer_heap.top().mds;
            }

            const auto before = cand_graph_set.size();
            {
                profile::ScopedPhase phase("extend_total");
                candidates_added_to_queue += extend(cand_graph_set, graph, cand_space, new_mds);
            }
            candidates_generated += static_cast<std::uint64_t>(cand_graph_set.size() - before);
        }

        const auto heap_threshold = static_cast<std::size_t>(static_cast<double>(k) * 1.5);
        while (cand_graph_set.size() < heap_threshold && !delayed_queue.empty()) {
            const auto tau_tth = get_tth_largest_ub(cand_graph_set, tth);
            const auto peek_mds = delayed_queue.peek_mds();
            if (!peek_mds.has_value()) {
                break;
            }
            if (*peek_mds > tau_tth) {
                const auto delayed = delayed_queue.pop_max();
                if (delayed.has_value()) {
                    profile::ScopedPhase phase("extend_total");
                    candidates_added_to_queue += extend_edge(cand_graph_set, delayed->second, delayed->first);
                }
            } else {
                delayed_queue.record_skipped();
                break;
            }
        }
    }

    answer_set = answer_heap.into_sorted_vec();
    std::reverse(answer_set.begin(), answer_set.end());
}

void MiningContext::print_answer() const {
    std::cout << "==================\n      Answer\n==================\n";
    std::cout << "Size of answer set: " << answer_set.size() << "\n";
    if (!answer_set.empty()) {
        std::cout << "Minimum MDS: " << answer_set.front().mds << "\n";
    }
    for (std::size_t i = 0; i < answer_set.size(); ++i) {
        std::cout << "--------\nPattern " << (i + 1) << ": MDS = " << answer_set[i].mds << "\n";
        for (const auto& code : answer_set[i].graph.dfs_code) {
            const auto fl_it = label_map_inverse.find(code.from_label);
            const auto tl_it = label_map_inverse.find(code.to_label);
            const auto el_it = edge_label_map_inverse.find(code.edge_label);
            const auto& fl = fl_it != label_map_inverse.end() ? fl_it->second : std::string{"?"};
            const auto& tl = tl_it != label_map_inverse.end() ? tl_it->second : std::string{"?"};
            const auto el = el_it != edge_label_map_inverse.end() ? el_it->second : -1;
            std::cout << "( " << code.from << ", " << code.to << ", " << fl << ", " << el << ", " << tl
                      << " )\n";
        }
    }
}

}
