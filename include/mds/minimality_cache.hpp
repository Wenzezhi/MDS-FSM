#pragma once

#include <cstdint>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "mds/dfs_code.hpp"
#include "mds/fxhash.hpp"

namespace mds {

class MinimalityCache {
public:
    explicit MinimalityCache(std::size_t max_size);

    static MinimalityCache with_default_size();
    static std::uint64_t hash_dfs_code(const std::vector<DfsCode>& dfs_code);

    template <typename Fn>
    bool check(std::uint64_t pattern_hash, Fn&& check_fn) {
        auto it = cache.find(pattern_hash);
        if (it != cache.end()) {
            ++hits;
            return it->second;
        }
        ++misses;
        const bool result = check_fn();
        if (cache.size() >= max_size_) {
            clear_half();
        }
        cache.emplace(pattern_hash, result);
        return result;
    }

    void clear_half();
    [[nodiscard]] double hit_rate() const;
    [[nodiscard]] std::size_t size() const { return cache.size(); }
    [[nodiscard]] std::size_t max_size() const { return max_size_; }
    void clear();
    [[nodiscard]] std::string stats() const;

    std::unordered_map<std::uint64_t, bool, fxhash::Hash<std::uint64_t>> cache;
    std::uint64_t hits = 0;
    std::uint64_t misses = 0;

private:
    std::size_t max_size_ = 0;
};

}
