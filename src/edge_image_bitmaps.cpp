#include "mds/edge_image_bitmaps.hpp"

#include <algorithm>
#include <vector>

#include "mds/simd_bitmap.hpp"

namespace mds {

EdgeImageBitmaps::EdgeImageBitmaps(std::size_t n_edges, std::size_t n_vertices)
    : bitmaps_(n_edges * 2, (n_vertices + 63) / 64),
      n_edges_(n_edges),
      n_vertices_(n_vertices) {}

void EdgeImageBitmaps::insert(std::size_t edge_id, std::size_t endpoint, std::size_t vertex) {
    bitmaps_.set_bit(edge_id * 2 + endpoint, vertex);
}

bool EdgeImageBitmaps::contains(std::size_t edge_id, std::size_t endpoint, std::size_t vertex) const {
    return bitmaps_.test_bit(edge_id * 2 + endpoint, vertex);
}

std::size_t EdgeImageBitmaps::size(std::size_t edge_id, std::size_t endpoint) const {
    return bitmaps_.popcount(edge_id * 2 + endpoint);
}

std::int32_t EdgeImageBitmaps::compute_mds(std::size_t edge_id) const {
    const auto row0 = edge_id * 2;
    const auto row1 = edge_id * 2 + 1;
    const auto img1_size = bitmaps_.popcount(row0);
    const auto img2_size = bitmaps_.popcount(row1);
    const auto union_size = compute_union_size(row0, row1);
    return static_cast<std::int32_t>(std::min({union_size / 2, img1_size, img2_size}));
}

std::size_t EdgeImageBitmaps::compute_union_size(std::size_t row0, std::size_t row1) const {
    const auto* lhs = bitmaps_.row(row0);
    const auto* rhs = bitmaps_.row(row1);
    const auto len = bitmaps_.row_capacity();
    std::size_t total = 0;
    for (std::size_t i = 0; i < len; ++i) {
#if defined(_MSC_VER)
        total += static_cast<std::size_t>(__popcnt64(lhs[i] | rhs[i]));
#else
        total += static_cast<std::size_t>(__builtin_popcountll(lhs[i] | rhs[i]));
#endif
    }
    return total;
}

}
