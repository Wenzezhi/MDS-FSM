#include "mds/simd_bitmap.hpp"

#include <algorithm>

#if defined(__AVX2__) || defined(_MSC_VER)
#include <immintrin.h>
#endif

namespace mds::simd_bitmap {

std::size_t popcount_scalar(const std::uint64_t* bitmap, std::size_t len) {
    std::size_t count = 0;
    for (std::size_t i = 0; i < len; ++i) {
#if defined(_MSC_VER)
        count += static_cast<std::size_t>(__popcnt64(bitmap[i]));
#else
        count += static_cast<std::size_t>(__builtin_popcountll(bitmap[i]));
#endif
    }
    return count;
}

void or_inplace_scalar(std::uint64_t* dst, const std::uint64_t* src, std::size_t len) {
    for (std::size_t i = 0; i < len; ++i) {
        dst[i] |= src[i];
    }
}

void and_inplace_scalar(std::uint64_t* dst, const std::uint64_t* src, std::size_t len) {
    for (std::size_t i = 0; i < len; ++i) {
        dst[i] &= src[i];
    }
}

void clear(std::uint64_t* bitmap, std::size_t len) {
    std::fill(bitmap, bitmap + len, 0ULL);
}

std::size_t popcount(const std::uint64_t* bitmap, std::size_t len) {
    return popcount_scalar(bitmap, len);
}

void or_inplace(std::uint64_t* dst, const std::uint64_t* src, std::size_t len) {
#if defined(__AVX2__)
    const std::size_t chunks = len / 4;
    for (std::size_t i = 0; i < chunks; ++i) {
        const auto idx = i * 4;
        const auto d = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(dst + idx));
        const auto s = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(src + idx));
        const auto result = _mm256_or_si256(d, s);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(dst + idx), result);
    }
    or_inplace_scalar(dst + chunks * 4, src + chunks * 4, len - chunks * 4);
#else
    or_inplace_scalar(dst, src, len);
#endif
}

void and_inplace(std::uint64_t* dst, const std::uint64_t* src, std::size_t len) {
#if defined(__AVX2__)
    const std::size_t chunks = len / 4;
    for (std::size_t i = 0; i < chunks; ++i) {
        const auto idx = i * 4;
        const auto d = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(dst + idx));
        const auto s = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(src + idx));
        const auto result = _mm256_and_si256(d, s);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(dst + idx), result);
    }
    and_inplace_scalar(dst + chunks * 4, src + chunks * 4, len - chunks * 4);
#else
    and_inplace_scalar(dst, src, len);
#endif
}

}
