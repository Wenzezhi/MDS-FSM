#include "mds/minimality_cache.hpp"

namespace mds {

MinimalityCache::MinimalityCache(std::size_t max_size) : max_size_(max_size) {
    cache.reserve(max_size_);
}

MinimalityCache MinimalityCache::with_default_size() {
    return MinimalityCache(262144);
}

std::uint64_t MinimalityCache::hash_dfs_code(const std::vector<DfsCode>& dfs_code) {
    fxhash::Hasher hasher;
    for (const auto& code : dfs_code) {
        fxhash::append(hasher, code.from);
        fxhash::append(hasher, code.to);
        fxhash::append(hasher, code.from_label);
        fxhash::append(hasher, code.edge_label);
        fxhash::append(hasher, code.to_label);
    }
    return static_cast<std::uint64_t>(hasher.finish());
}

void MinimalityCache::clear_half() {
    cache.clear();
    cache.reserve(max_size_);
}

double MinimalityCache::hit_rate() const {
    const auto total = hits + misses;
    if (total == 0) {
        return 0.0;
    }
    return static_cast<double>(hits) / static_cast<double>(total) * 100.0;
}

void MinimalityCache::clear() {
    cache.clear();
    cache.reserve(max_size_);
    hits = 0;
    misses = 0;
}

std::string MinimalityCache::stats() const {
    std::ostringstream oss;
    oss << "MinimalityCache: " << size() << " entries, " << hit_rate() << "% hit rate ("
        << hits << " hits, " << misses << " misses)";
    return oss.str();
}

}
