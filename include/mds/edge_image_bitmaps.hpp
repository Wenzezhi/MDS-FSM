#pragma once

#include <cstddef>
#include <cstdint>

#include "mds/flat_bitmap.hpp"

namespace mds {

class EdgeImageBitmaps {
public:
    EdgeImageBitmaps() = default;
    EdgeImageBitmaps(std::size_t n_edges, std::size_t n_vertices);

    void insert(std::size_t edge_id, std::size_t endpoint, std::size_t vertex);
    [[nodiscard]] bool contains(std::size_t edge_id, std::size_t endpoint, std::size_t vertex) const;
    [[nodiscard]] std::size_t size(std::size_t edge_id, std::size_t endpoint) const;
    [[nodiscard]] std::int32_t compute_mds(std::size_t edge_id) const;
    [[nodiscard]] std::size_t n_edges() const { return n_edges_; }
    [[nodiscard]] std::size_t n_vertices() const { return n_vertices_; }

private:
    [[nodiscard]] std::size_t compute_union_size(std::size_t row0, std::size_t row1) const;

    FlatBitmap bitmaps_;
    std::size_t n_edges_ = 0;
    std::size_t n_vertices_ = 0;
};

}
