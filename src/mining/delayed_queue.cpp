#include "mds/mining/delayed_queue.hpp"

namespace mds {

DelayedExtensionQueue::DelayedExtensionQueue(std::size_t tth) : tth_(tth) {}

void DelayedExtensionQueue::push(std::int32_t mds, Graph pattern) {
    queue_.push(DelayedElement{mds, std::move(pattern)});
    ++delayed_count_;
}

std::optional<std::pair<std::int32_t, Graph>> DelayedExtensionQueue::pop_max() {
    if (queue_.empty()) {
        return std::nullopt;
    }
    auto elem = std::move(const_cast<DelayedElement&>(queue_.top()));
    queue_.pop();
    return std::make_pair(elem.mds, std::move(elem.pattern));
}

std::optional<std::int32_t> DelayedExtensionQueue::peek_mds() const {
    if (queue_.empty()) {
        return std::nullopt;
    }
    return queue_.top().mds;
}

void DelayedExtensionQueue::record_immediate() {
    ++immediate_count_;
}

void DelayedExtensionQueue::record_skipped() {
    ++skipped_count_;
}

std::tuple<std::size_t, std::size_t, std::size_t> DelayedExtensionQueue::stats() const {
    return {immediate_count_, delayed_count_, skipped_count_};
}

}
