#pragma once

#include <algorithm>
#include <cstdint>
#include <memory>
#include <vector>

#include "mds/candidate_space.hpp"
#include "mds/graph.hpp"

namespace mds {

struct HeapPayload {
    Graph graph;
    std::vector<CandidateSpace> cand_space;
};

struct HeapElement {
    std::int32_t ub = 0;
    std::int32_t nv = 0;
    bool is_backward = false;
    std::unique_ptr<HeapPayload> payload;

    HeapElement() = default;
    HeapElement(HeapElement&&) noexcept = default;
    HeapElement& operator=(HeapElement&&) noexcept = default;
    HeapElement(const HeapElement&) = delete;
    HeapElement& operator=(const HeapElement&) = delete;

    [[nodiscard]] static HeapElement from_payload(std::int32_t ub,
                                                  std::int32_t nv,
                                                  bool is_backward,
                                                  Graph&& graph,
                                                  std::vector<CandidateSpace>&& cand_space) {
        HeapElement element;
        element.ub = ub;
        element.nv = nv;
        element.is_backward = is_backward;
        element.payload = std::make_unique<HeapPayload>(HeapPayload{
            .graph = std::move(graph),
            .cand_space = std::move(cand_space),
        });
        return element;
    }
};

struct HeapElementLess {
    bool operator()(const HeapElement& lhs, const HeapElement& rhs) const {
        if (lhs.ub != rhs.ub) {
            return lhs.ub < rhs.ub;
        }
        if (lhs.nv != rhs.nv) {
            return rhs.nv < lhs.nv;
        }
        return lhs.is_backward < rhs.is_backward;
    }
};

class CandidateHeap {
public:
    [[nodiscard]] bool empty() const noexcept { return data_.empty(); }
    [[nodiscard]] std::size_t size() const noexcept { return data_.size(); }
    [[nodiscard]] const HeapElement& top() const { return data_.front(); }
    [[nodiscard]] const std::vector<HeapElement>& data() const noexcept { return data_; }

    void push(HeapElement&& item);
    void pop();
    [[nodiscard]] HeapElement pop_value();

private:
    std::vector<HeapElement> data_;
};

std::int32_t get_tth_largest_ub(const CandidateHeap& heap, std::size_t tth);
HeapElement pop_heap_element(CandidateHeap& heap);

struct AnswerElement {
    std::int32_t mds = 0;
    Graph graph;
};

class AnswerHeap {
public:
    [[nodiscard]] bool empty() const noexcept { return data_.empty(); }
    [[nodiscard]] std::size_t size() const noexcept { return data_.size(); }
    [[nodiscard]] const AnswerElement& top() const { return data_.front(); }

    void push(const AnswerElement& item);
    void push(AnswerElement&& item);
    void pop();
    [[nodiscard]] AnswerElement pop_value();
    [[nodiscard]] std::vector<AnswerElement> into_sorted_vec();

private:
    static std::size_t sift_up(std::vector<AnswerElement>& data, std::size_t start, std::size_t pos);
    static std::size_t sift_down_range(std::vector<AnswerElement>& data, std::size_t pos, std::size_t end);
    static void sift_down_to_bottom(std::vector<AnswerElement>& data, std::size_t pos);

    std::vector<AnswerElement> data_;
};
AnswerElement pop_answer_element(AnswerHeap& heap);

}
