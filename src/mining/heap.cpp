#include "mds/mining/heap.hpp"

namespace mds {

std::int32_t get_tth_largest_ub(const CandidateHeap& heap, std::size_t tth) {
    if (heap.empty()) {
        return 0;
    }

    if (tth >= heap.size()) {
        auto min_value = heap.top().ub;
        for (const auto& element : heap.data()) {
            min_value = std::min(min_value, element.ub);
        }
        return min_value;
    }

    std::vector<std::int32_t> values;
    values.reserve(heap.size());
    for (const auto& element : heap.data()) {
        values.push_back(element.ub);
    }
    std::sort(values.begin(), values.end(), std::greater<>{});
    return values[tth - 1];
}

HeapElement pop_heap_element(CandidateHeap& heap) {
    return heap.pop_value();
}

void CandidateHeap::push(HeapElement&& item) {
    data_.push_back(std::move(item));
    std::push_heap(data_.begin(), data_.end(), HeapElementLess{});
}

void CandidateHeap::pop() {
    (void)pop_value();
}

HeapElement CandidateHeap::pop_value() {
    std::pop_heap(data_.begin(), data_.end(), HeapElementLess{});
    auto item = std::move(data_.back());
    data_.pop_back();
    return item;
}

AnswerElement pop_answer_element(AnswerHeap& heap) {
    return heap.pop_value();
}

std::size_t AnswerHeap::sift_up(std::vector<AnswerElement>& data, std::size_t start, std::size_t pos) {
    auto item = std::move(data[pos]);
    while (pos > start) {
        const auto parent = (pos - 1) / 2;
        if (item.mds >= data[parent].mds) {
            break;
        }
        data[pos] = std::move(data[parent]);
        pos = parent;
    }
    data[pos] = std::move(item);
    return pos;
}

std::size_t AnswerHeap::sift_down_range(std::vector<AnswerElement>& data, std::size_t pos, std::size_t end) {
    auto item = std::move(data[pos]);
    auto child = 2 * pos + 1;

    while (child + 1 < end) {
        child += (data[child].mds >= data[child + 1].mds) ? 1U : 0U;
        if (item.mds <= data[child].mds) {
            data[pos] = std::move(item);
            return pos;
        }
        data[pos] = std::move(data[child]);
        pos = child;
        child = 2 * pos + 1;
    }

    if (child == end - 1 && item.mds > data[child].mds) {
        data[pos] = std::move(data[child]);
        pos = child;
    }

    data[pos] = std::move(item);
    return pos;
}

void AnswerHeap::sift_down_to_bottom(std::vector<AnswerElement>& data, std::size_t pos) {
    const auto end = data.size();
    const auto start = pos;
    auto item = std::move(data[pos]);
    auto child = 2 * pos + 1;

    while (child + 1 < end) {
        child += (data[child].mds >= data[child + 1].mds) ? 1U : 0U;
        data[pos] = std::move(data[child]);
        pos = child;
        child = 2 * pos + 1;
    }

    if (child == end - 1) {
        data[pos] = std::move(data[child]);
        pos = child;
    }

    data[pos] = std::move(item);
    sift_up(data, start, pos);
}

void AnswerHeap::push(const AnswerElement& item) {
    const auto old_len = data_.size();
    data_.push_back(item);
    sift_up(data_, 0, old_len);
}

void AnswerHeap::push(AnswerElement&& item) {
    const auto old_len = data_.size();
    data_.push_back(std::move(item));
    sift_up(data_, 0, old_len);
}

void AnswerHeap::pop() {
    (void)pop_value();
}

AnswerElement AnswerHeap::pop_value() {
    auto item = std::move(data_.back());
    data_.pop_back();
    if (!data_.empty()) {
        std::swap(item, data_[0]);
        sift_down_to_bottom(data_, 0);
    }
    return item;
}

std::vector<AnswerElement> AnswerHeap::into_sorted_vec() {
    auto data = std::move(data_);
    auto end = data.size();
    while (end > 1) {
        --end;
        std::swap(data[0], data[end]);
        sift_down_range(data, 0, end);
    }
    return data;
}

}
