#pragma once

#include <cstddef>
#include <cstdint>

namespace mds {

using NumberType = std::int32_t;

inline constexpr std::size_t kMaxNumVertex = 257;
inline constexpr std::size_t kBitsPerLabel = 4;
inline constexpr std::size_t kUint64Size = 64;
inline constexpr std::size_t kMaxQueryVertex = 32;
inline constexpr std::size_t kSmallQueryDegree = 32;

}
