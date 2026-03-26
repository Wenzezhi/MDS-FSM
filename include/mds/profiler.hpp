#pragma once

#include <chrono>
#include <cstdint>
#include <string_view>

namespace mds::profile {

struct PhaseSnapshot {
    std::uint64_t calls = 0;
    std::uint64_t total_ns = 0;
    std::uint64_t max_ns = 0;
};

bool is_enabled();
void add_duration(std::string_view phase_name, std::uint64_t elapsed_ns);
void set_counter(std::string_view counter_name, std::uint64_t value);
void write_to_configured_path();

class ScopedPhase {
public:
    explicit ScopedPhase(std::string_view phase_name);
    ~ScopedPhase();

    ScopedPhase(const ScopedPhase&) = delete;
    ScopedPhase& operator=(const ScopedPhase&) = delete;

private:
    std::string_view phase_name_{};
    std::chrono::steady_clock::time_point start_{};
    bool active_ = false;
};

}
