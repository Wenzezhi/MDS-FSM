#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <queue>
#include <utility>
#include <vector>

#include "mds/graph.hpp"

namespace mds {

struct DelayedElement {
    std::int32_t mds = 0;
    Graph pattern;
};

struct DelayedElementLess {
    bool operator()(const DelayedElement& lhs, const DelayedElement& rhs) const {
        return lhs.mds < rhs.mds;
    }
};

class DelayedExtensionQueue {
public:
    explicit DelayedExtensionQueue(std::size_t tth);

    void push(std::int32_t mds, Graph pattern);
    [[nodiscard]] std::optional<std::pair<std::int32_t, Graph>> pop_max();
    [[nodiscard]] std::optional<std::int32_t> peek_mds() const;
    [[nodiscard]] bool empty() const { return queue_.empty(); }
    [[nodiscard]] std::size_t len() const { return queue_.size(); }
    [[nodiscard]] std::size_t get_tth() const { return tth_; }
    void record_immediate();
    void record_skipped();
    [[nodiscard]] std::tuple<std::size_t, std::size_t, std::size_t> stats() const;

private:
    std::priority_queue<DelayedElement, std::vector<DelayedElement>, DelayedElementLess> queue_;
    std::size_t tth_ = 0;
    std::size_t immediate_count_ = 0;
    std::size_t delayed_count_ = 0;
    std::size_t skipped_count_ = 0;
};

}
