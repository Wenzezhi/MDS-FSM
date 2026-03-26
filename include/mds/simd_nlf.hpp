#pragma once

#include <cstddef>
#include <cstdint>

namespace mds {

class SimdNlfComparator {
public:
    explicit SimdNlfComparator(std::size_t nlf_size) : nlf_size_(nlf_size) {}

    [[nodiscard]] bool check(const std::uint64_t* pattern_nlf, const std::uint64_t* data_nlf) const;
    [[nodiscard]] bool check_scalar(const std::uint64_t* pattern_nlf, const std::uint64_t* data_nlf) const;

private:
    std::size_t nlf_size_ = 0;
};

}
