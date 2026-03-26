#pragma once

#include <algorithm>
#include <cstddef>
#include <vector>

namespace mds {

template <typename T>
class FlatMatrix {
public:
    FlatMatrix() = default;
    FlatMatrix(std::size_t num_rows, std::size_t row_capacity)
        : data_(num_rows * row_capacity, T{}),
          row_capacity_(row_capacity),
          num_rows_(num_rows) {}

    [[nodiscard]] T get(std::size_t row, std::size_t col) const {
        return data_[row * row_capacity_ + col];
    }

    [[nodiscard]] const T& at(std::size_t row, std::size_t col) const {
        return data_[row * row_capacity_ + col];
    }

    void set(std::size_t row, std::size_t col, T value) {
        data_[row * row_capacity_ + col] = value;
    }

    [[nodiscard]] T& at_mut(std::size_t row, std::size_t col) {
        return data_[row * row_capacity_ + col];
    }

    [[nodiscard]] T* get_mut(std::size_t row, std::size_t col) {
        return &data_[row * row_capacity_ + col];
    }

    [[nodiscard]] T* row_mut(std::size_t row) {
        return data_.data() + row * row_capacity_;
    }

    [[nodiscard]] const T* row(std::size_t row) const {
        return data_.data() + row * row_capacity_;
    }

    [[nodiscard]] std::size_t num_rows() const { return num_rows_; }
    [[nodiscard]] std::size_t row_capacity() const { return row_capacity_; }
    [[nodiscard]] std::size_t len() const { return num_rows_; }
    [[nodiscard]] bool empty() const { return num_rows_ == 0; }

    void fill(const T& value) {
        std::fill(data_.begin(), data_.end(), value);
    }

    void swap_rows(std::size_t row1, std::size_t row2) {
        if (row1 == row2) {
            return;
        }
        const auto start1 = row1 * row_capacity_;
        const auto start2 = row2 * row_capacity_;
        for (std::size_t i = 0; i < row_capacity_; ++i) {
            std::swap(data_[start1 + i], data_[start2 + i]);
        }
    }

    void copy_row(std::size_t dst_row, std::size_t src_row) {
        if (dst_row == src_row) {
            return;
        }
        std::copy_n(row(src_row), row_capacity_, row_mut(dst_row));
    }

private:
    std::vector<T> data_;
    std::size_t row_capacity_ = 0;
    std::size_t num_rows_ = 0;
};

}
