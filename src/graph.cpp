#include "mds/graph.hpp"

#include <algorithm>

namespace mds {

std::size_t Graph::pair_matrix_index(std::size_t u, std::size_t v) const {
    return u * static_cast<std::size_t>(n_vertex) + v;
}

void Graph::init_vertex(NumberType n) {
    n_vertex = n;
    const auto nu = static_cast<std::size_t>(n);
    label.assign(nu, 0);
    nbr_offset.assign(nu + 1, 0);
    core.assign(nu, 0);
    degree.assign(nu, 0);
    vertex.assign(nu, 0);
    max_nbr_degree.assign(nu, 0);
    right_most_path.assign(nu, 0);
}

void Graph::init_nbr_to_pos() {
    const auto nu = static_cast<std::size_t>(n_vertex);
    nbr_to_pos_data.assign(nu * nu, -1);
    for (std::size_t i = 0; i < nu; ++i) {
        for (auto j = static_cast<std::size_t>(nbr_offset[i]); j < static_cast<std::size_t>(nbr_offset[i + 1]); ++j) {
            const auto un = static_cast<std::size_t>(nbr[j].first);
            set_nbr_to_pos(i, un, static_cast<std::int32_t>(j - static_cast<std::size_t>(nbr_offset[i])));
        }
    }
}

void Graph::init_edge(NumberType n) {
    n_edge = n;
    nbr.assign(static_cast<std::size_t>(n) * 2, {0, 0});
}

void Graph::sort_neighbors(std::size_t n_unique_label) {
    vertex_offset.assign(n_unique_label + 1, 0);
    label_nbr_offset_stride = n_unique_label;
    label_nbr_offset_table.assign(static_cast<std::size_t>(n_vertex) * n_unique_label, {0, 0});

    for (std::size_t i = 0; i < static_cast<std::size_t>(n_vertex); ++i) {
        ++vertex_offset[static_cast<std::size_t>(label[i])];
        if (degree[i] == 0) {
            continue;
        }

        const auto start = static_cast<std::size_t>(nbr_offset[i]);
        const auto end = static_cast<std::size_t>(nbr_offset[i + 1]);
        std::sort(nbr.begin() + static_cast<std::ptrdiff_t>(start),
                  nbr.begin() + static_cast<std::ptrdiff_t>(end),
                  [&](const auto& a, const auto& b) {
                      const auto la = label[static_cast<std::size_t>(a.first)];
                      const auto lb = label[static_cast<std::size_t>(b.first)];
                      if (la != lb) {
                          return la < lb;
                      }
                      const auto da = degree[static_cast<std::size_t>(a.first)];
                      const auto db = degree[static_cast<std::size_t>(b.first)];
                      if (da != db) {
                          return da < db;
                      }
                      return a.first < b.first;
                  });

        NumberType curr_cnt = 1;
        NumberType curr_index = static_cast<NumberType>(start);
        auto prev = label[static_cast<std::size_t>(nbr[start].first)];

        if (end - start > 1) {
            for (auto j = start + 1; j < end; ++j) {
                const auto curr = label[static_cast<std::size_t>(nbr[j].first)];
                if (prev == curr) {
                    ++curr_cnt;
                } else {
                    label_nbr_offset_table[i * label_nbr_offset_stride + static_cast<std::size_t>(prev)] = {curr_index,
                                                                                                             curr_cnt};
                    curr_index += curr_cnt;
                    prev = curr;
                    curr_cnt = 1;
                }
            }
        }
        label_nbr_offset_table[i * label_nbr_offset_stride + static_cast<std::size_t>(prev)] = {curr_index, curr_cnt};
    }

    auto curr = vertex_offset[0];
    vertex_offset[0] = 0;
    for (std::size_t i = 1; i <= n_unique_label; ++i) {
        const auto prev = vertex_offset[i];
        vertex_offset[i] = vertex_offset[i - 1] + curr;
        curr = prev;
    }

    for (std::size_t i = 0; i < static_cast<std::size_t>(n_vertex); ++i) {
        const auto l = static_cast<std::size_t>(label[i]);
        vertex[static_cast<std::size_t>(vertex_offset[l])] = static_cast<std::int32_t>(i);
        ++vertex_offset[l];
    }

    for (std::size_t i = n_unique_label; i > 0; --i) {
        vertex_offset[i] = vertex_offset[i - 1];
    }
    vertex_offset[0] = 0;

    for (std::size_t i = 0; i < n_unique_label; ++i) {
        const auto start = static_cast<std::size_t>(vertex_offset[i]);
        const auto end = static_cast<std::size_t>(vertex_offset[i + 1]);
        std::sort(vertex.begin() + static_cast<std::ptrdiff_t>(start),
                  vertex.begin() + static_cast<std::ptrdiff_t>(end),
                  [&](std::int32_t a, std::int32_t b) {
                      const auto da = degree[static_cast<std::size_t>(a)];
                      const auto db = degree[static_cast<std::size_t>(b)];
                      if (da != db) {
                          return da < db;
                      }
                      return a < b;
                  });
    }
}

void Graph::compute_nlf(std::size_t nlf_size, std::size_t n_unique_label) {
    nlf.assign(static_cast<std::size_t>(n_vertex) * nlf_size, 0ULL);
    std::vector<std::int32_t> cnt_label(n_unique_label, 0);
    std::vector<std::size_t> touched_labels(n_unique_label, 0);

    for (std::size_t j = 0; j < static_cast<std::size_t>(n_vertex); ++j) {
        if (degree[j] == 0) {
            continue;
        }
        std::size_t touched_count = 0;

        for (auto k = static_cast<std::size_t>(nbr_offset[j]); k < static_cast<std::size_t>(nbr_offset[j + 1]); ++k) {
            const auto neighbor = static_cast<std::size_t>(nbr[k].first);
            if (degree[neighbor] > max_nbr_degree[j]) {
                max_nbr_degree[j] = degree[neighbor];
            }
            const auto lid = static_cast<std::size_t>(label[neighbor]);
            if (cnt_label[lid] == 0) {
                touched_labels[touched_count++] = lid;
            }
            if (static_cast<std::size_t>(cnt_label[lid]) < kBitsPerLabel) {
                const auto start_bit = lid * kBitsPerLabel;
                const auto idx = nlf_size - 1 - (start_bit + static_cast<std::size_t>(cnt_label[lid])) / kUint64Size;
                const auto pos = (start_bit + static_cast<std::size_t>(cnt_label[lid])) % kUint64Size;
                nlf[j * nlf_size + idx] |= (1ULL << pos);
                ++cnt_label[lid];
            }
        }

        for (std::size_t idx = 0; idx < touched_count; ++idx) {
            cnt_label[touched_labels[idx]] = 0;
        }
    }
}

void Graph::compute_nbr_safety() {
    comp_label_idx_data.assign(static_cast<std::size_t>(n_vertex) * static_cast<std::size_t>(n_vertex), 0);
    nbr_label_data.assign(nbr.size(), 0);
    nbr_comp_idx_data.assign(nbr.size(), 0);
    nbr_reverse_pos_data.assign(nbr.size(), -1);
    nbr_reverse_comp_idx_data.assign(nbr.size(), 0);

    for (std::size_t u = 0; u < static_cast<std::size_t>(n_vertex); ++u) {
        std::unordered_map<std::pair<std::int32_t, std::int32_t>,
                           std::int32_t,
                           fxhash::Hash<std::pair<std::int32_t, std::int32_t>>> vis;
        for (auto j = static_cast<std::size_t>(nbr_offset[u]); j < static_cast<std::size_t>(nbr_offset[u + 1]); ++j) {
            const auto un = static_cast<std::size_t>(nbr[j].first);
            const auto qel = nbr[j].second;
            const auto lun = label[un];
            const auto key = std::make_pair(lun, qel);
            auto it = vis.find(key);
            std::int32_t idx = 0;
            if (it != vis.end()) {
                idx = it->second;
            } else {
                idx = static_cast<std::int32_t>(vis.size());
                vis.emplace(key, idx);
            }
            set_comp_label_idx(u, un, idx);
            nbr_label_data[j] = lun;
            nbr_comp_idx_data[j] = idx;
        }
    }

    for (std::size_t u = 0; u < static_cast<std::size_t>(n_vertex); ++u) {
        for (auto j = static_cast<std::size_t>(nbr_offset[u]); j < static_cast<std::size_t>(nbr_offset[u + 1]); ++j) {
            const auto un = static_cast<std::size_t>(nbr[j].first);
            nbr_reverse_pos_data[j] = nbr_to_pos_at(un, u);
            nbr_reverse_comp_idx_data[j] = comp_label_idx_at(un, u);
        }
    }
}

std::pair<NumberType, NumberType> Graph::lookup_label_nbr_offset(std::int32_t vertex_id, std::int32_t label_id) const {
    if (vertex_id < 0 || label_id < 0) {
        return {0, 0};
    }
    const auto vertex = static_cast<std::size_t>(vertex_id);
    const auto label_idx = static_cast<std::size_t>(label_id);
    if (vertex >= static_cast<std::size_t>(n_vertex) || label_idx >= label_nbr_offset_stride) {
        return {0, 0};
    }
    return label_nbr_offset_table[vertex * label_nbr_offset_stride + label_idx];
}

std::int32_t Graph::nbr_to_pos_at(std::size_t u, std::size_t v) const {
    return nbr_to_pos_data[pair_matrix_index(u, v)];
}

void Graph::set_nbr_to_pos(std::size_t u, std::size_t v, std::int32_t value) {
    nbr_to_pos_data[pair_matrix_index(u, v)] = value;
}

std::int32_t Graph::comp_label_idx_at(std::size_t u, std::size_t v) const {
    return comp_label_idx_data[pair_matrix_index(u, v)];
}

void Graph::set_comp_label_idx(std::size_t u, std::size_t v, std::int32_t value) {
    comp_label_idx_data[pair_matrix_index(u, v)] = value;
}

std::uint64_t hash_pattern(const std::vector<DfsCode>& dfs_code) {
    fxhash::Hasher hasher;
    fxhash::append(hasher, dfs_code);
    return static_cast<std::uint64_t>(hasher.finish());
}

}
