#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace mds::fxhash {

class Hasher {
public:
    Hasher() = default;

    void add_u64(std::uint64_t value) {
        hash_ = rotl(hash_, 5) ^ value;
        hash_ *= kConst;
    }

    void add_bytes(std::string_view bytes) {
        std::size_t i = 0;
        const auto* ptr = reinterpret_cast<const unsigned char*>(bytes.data());
        while (i + sizeof(std::uint64_t) <= bytes.size()) {
            std::uint64_t word = 0;
            for (std::size_t j = 0; j < sizeof(std::uint64_t); ++j) {
                word |= static_cast<std::uint64_t>(ptr[i + j]) << (j * 8);
            }
            add_u64(word);
            i += sizeof(std::uint64_t);
        }
        if (i + sizeof(std::uint32_t) <= bytes.size()) {
            std::uint32_t word = 0;
            for (std::size_t j = 0; j < sizeof(std::uint32_t); ++j) {
                word |= static_cast<std::uint32_t>(ptr[i + j]) << (j * 8);
            }
            add_u64(word);
            i += sizeof(std::uint32_t);
        }
        if (i + sizeof(std::uint16_t) <= bytes.size()) {
            std::uint16_t word = static_cast<std::uint16_t>(ptr[i]) |
                                 static_cast<std::uint16_t>(ptr[i + 1] << 8);
            add_u64(word);
            i += sizeof(std::uint16_t);
        }
        if (i < bytes.size()) {
            add_u64(ptr[i]);
        }
    }

    [[nodiscard]] std::size_t finish() const {
        return static_cast<std::size_t>(hash_);
    }

private:
    static constexpr std::uint64_t kConst = 0x517cc1b727220a95ULL;
    std::uint64_t hash_ = 0;

    static constexpr std::uint64_t rotl(std::uint64_t value, int shift) {
        return (value << shift) | (value >> (64 - shift));
    }
};

template <typename T>
inline void append(Hasher& hasher, const T& value);

template <typename T, bool IsEnum = std::is_enum_v<T>>
struct UnsignedLike;

template <typename T>
struct UnsignedLike<T, false> {
    using type = std::make_unsigned_t<T>;
};

template <typename T>
struct UnsignedLike<T, true> {
    using type = std::make_unsigned_t<std::underlying_type_t<T>>;
};

template <typename T>
struct Hash {
    std::size_t operator()(const T& value) const {
        Hasher hasher;
        append(hasher, value);
        return hasher.finish();
    }
};

template <typename T>
requires std::is_integral_v<T> || std::is_enum_v<T>
inline void append(Hasher& hasher, const T& value) {
    using Unsigned = typename UnsignedLike<T>::type;
    hasher.add_u64(static_cast<std::uint64_t>(static_cast<Unsigned>(value)));
}

inline void append(Hasher& hasher, const std::string& value) {
    hasher.add_bytes(value);
}

inline void append(Hasher& hasher, std::string_view value) {
    hasher.add_bytes(value);
}

template <typename A, typename B>
inline void append(Hasher& hasher, const std::pair<A, B>& value) {
    append(hasher, value.first);
    append(hasher, value.second);
}

template <typename... Ts, std::size_t... Is>
inline void append_tuple_impl(Hasher& hasher,
                              const std::tuple<Ts...>& value,
                              std::index_sequence<Is...>) {
    (append(hasher, std::get<Is>(value)), ...);
}

template <typename... Ts>
inline void append(Hasher& hasher, const std::tuple<Ts...>& value) {
    append_tuple_impl(hasher, value, std::index_sequence_for<Ts...>{});
}

template <typename T>
inline void append(Hasher& hasher, const std::vector<T>& values) {
    for (const auto& value : values) {
        append(hasher, value);
    }
}

}
