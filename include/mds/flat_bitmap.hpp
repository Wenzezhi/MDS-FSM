#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace mds {

class FlatBitmap {
public:
    FlatBitmap() = default;
    FlatBitmap(std::size_t num_rows, std::size_t row_capacity);

    [[nodiscard]] const std::uint64_t* row(std::size_t row_index) const;
    [[nodiscard]] std::uint64_t* row_mut(std::size_t row_index);

    [[nodiscard]] std::size_t num_rows() const { return num_rows_; }
    [[nodiscard]] std::size_t row_capacity() const { return row_capacity_; }

    void set_bit(std::size_t row_index, std::size_t bit_index);
    void clear_bit(std::size_t row_index, std::size_t bit_index);
    [[nodiscard]] bool test_bit(std::size_t row_index, std::size_t bit_index) const;
    void clear_row(std::size_t row_index);
    [[nodiscard]] std::size_t popcount(std::size_t row_index) const;
    void or_inplace(std::size_t dst_row, std::size_t src_row);
    void and_inplace(std::size_t dst_row, std::size_t src_row);

private:
    std::vector<std::uint64_t> data_;
    std::size_t row_capacity_ = 0;
    std::size_t num_rows_ = 0;
};

}
