#include "mds/simd_nlf.hpp"

#if defined(__AVX2__) || defined(_MSC_VER)
#include <immintrin.h>
#endif

namespace mds {

bool SimdNlfComparator::check_scalar(const std::uint64_t* pattern_nlf,
                                     const std::uint64_t* data_nlf) const {
    for (std::size_t i = 0; i < nlf_size_; ++i) {
        if (data_nlf[i] != (pattern_nlf[i] | data_nlf[i])) {
            return false;
        }
    }
    return true;
}

bool SimdNlfComparator::check(const std::uint64_t* pattern_nlf, const std::uint64_t* data_nlf) const {
#if defined(__AVX2__)
    std::size_t i = 0;
    while (i + 4 <= nlf_size_) {
        const auto pattern = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(pattern_nlf + i));
        const auto data = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(data_nlf + i));
        const auto or_result = _mm256_or_si256(pattern, data);
        const auto cmp = _mm256_cmpeq_epi64(data, or_result);
        const auto mask = _mm256_movemask_pd(_mm256_castsi256_pd(cmp));
        if (mask != 0xF) {
            return false;
        }
        i += 4;
    }
    for (; i < nlf_size_; ++i) {
        if (data_nlf[i] != (pattern_nlf[i] | data_nlf[i])) {
            return false;
        }
    }
    return true;
#else
    return check_scalar(pattern_nlf, data_nlf);
#endif
}

}
