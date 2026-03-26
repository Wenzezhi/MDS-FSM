#include "mds/flat_bitmap.hpp"

#include "mds/simd_bitmap.hpp"

namespace mds {

FlatBitmap::FlatBitmap(std::size_t num_rows, std::size_t row_capacity)
    : data_(num_rows * row_capacity, 0ULL),
      row_capacity_(row_capacity),
      num_rows_(num_rows) {}

const std::uint64_t* FlatBitmap::row(std::size_t row_index) const {
    return data_.data() + row_index * row_capacity_;
}

std::uint64_t* FlatBitmap::row_mut(std::size_t row_index) {
    return data_.data() + row_index * row_capacity_;
}

void FlatBitmap::set_bit(std::size_t row_index, std::size_t bit_index) {
    const auto offset = row_index * row_capacity_ + bit_index / 64;
    const auto bit = bit_index % 64;
    data_[offset] |= (1ULL << bit);
}

void FlatBitmap::clear_bit(std::size_t row_index, std::size_t bit_index) {
    const auto offset = row_index * row_capacity_ + bit_index / 64;
    const auto bit = bit_index % 64;
    data_[offset] &= ~(1ULL << bit);
}

bool FlatBitmap::test_bit(std::size_t row_index, std::size_t bit_index) const {
    const auto offset = row_index * row_capacity_ + bit_index / 64;
    const auto bit = bit_index % 64;
    return ((data_[offset] >> bit) & 1ULL) == 1ULL;
}

void FlatBitmap::clear_row(std::size_t row_index) {
    simd_bitmap::clear(row_mut(row_index), row_capacity_);
}

std::size_t FlatBitmap::popcount(std::size_t row_index) const {
    return simd_bitmap::popcount(row(row_index), row_capacity_);
}

void FlatBitmap::or_inplace(std::size_t dst_row, std::size_t src_row) {
    if (dst_row == src_row) {
        return;
    }
    simd_bitmap::or_inplace(row_mut(dst_row), row(src_row), row_capacity_);
}

void FlatBitmap::and_inplace(std::size_t dst_row, std::size_t src_row) {
    if (dst_row == src_row) {
        return;
    }
    simd_bitmap::and_inplace(row_mut(dst_row), row(src_row), row_capacity_);
}

}
