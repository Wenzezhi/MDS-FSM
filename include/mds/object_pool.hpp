#pragma once

#include <unordered_set>
#include <vector>

namespace mds {

template <typename T>
class ObjectPool {
public:
    explicit ObjectPool(std::size_t capacity) : in_use_(0) {
        available_.reserve(capacity);
        for (std::size_t i = 0; i < capacity; ++i) {
            available_.push_back(T{});
        }
    }

    [[nodiscard]] T acquire() {
        ++in_use_;
        if (available_.empty()) {
            return T{};
        }
        auto value = std::move(available_.back());
        available_.pop_back();
        return value;
    }

    void release(T value) {
        --in_use_;
        available_.push_back(std::move(value));
    }

    void reset() const {}

private:
    std::vector<T> available_;
    std::size_t in_use_;
};

template <typename T>
inline void release_vec(ObjectPool<std::vector<T>>& pool, std::vector<T> value) {
    value.clear();
    pool.release(std::move(value));
}

template <typename T, typename Hash, typename Eq>
inline void release_set(ObjectPool<std::unordered_set<T, Hash, Eq>>& pool,
                        std::unordered_set<T, Hash, Eq> value) {
    value.clear();
    pool.release(std::move(value));
}

}
