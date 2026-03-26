#pragma once

#include <cstddef>
#include <cstdint>

namespace mds::simd_bitmap {

std::size_t popcount_scalar(const std::uint64_t* bitmap, std::size_t len);
void or_inplace_scalar(std::uint64_t* dst, const std::uint64_t* src, std::size_t len);
void and_inplace_scalar(std::uint64_t* dst, const std::uint64_t* src, std::size_t len);
void clear(std::uint64_t* bitmap, std::size_t len);
std::size_t popcount(const std::uint64_t* bitmap, std::size_t len);
void or_inplace(std::uint64_t* dst, const std::uint64_t* src, std::size_t len);
void and_inplace(std::uint64_t* dst, const std::uint64_t* src, std::size_t len);

}
