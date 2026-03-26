#include "mds/profiler.hpp"

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <mutex>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace mds::profile {

namespace {

struct ProfilerState {
    bool initialized = false;
    bool enabled = false;
    std::string output_path;
    std::unordered_map<std::string, PhaseSnapshot> phases;
    std::unordered_map<std::string, std::uint64_t> counters;
    std::mutex mutex;
};

ProfilerState& state() {
    static ProfilerState profiler_state;
    if (!profiler_state.initialized) {
        profiler_state.initialized = true;
        if (const auto* value = std::getenv("MDS_PROFILE_OUT")) {
            profiler_state.output_path = value;
            profiler_state.enabled = !profiler_state.output_path.empty();
        }
    }
    return profiler_state;
}

std::string escape_json(std::string_view value) {
    std::string escaped;
    escaped.reserve(value.size() + 8);
    for (const auto ch : value) {
        switch (ch) {
            case '\\':
                escaped += "\\\\";
                break;
            case '"':
                escaped += "\\\"";
                break;
            case '\n':
                escaped += "\\n";
                break;
            case '\r':
                escaped += "\\r";
                break;
            case '\t':
                escaped += "\\t";
                break;
            default:
                escaped += ch;
                break;
        }
    }
    return escaped;
}

}

bool is_enabled() {
    return state().enabled;
}

void add_duration(std::string_view phase_name, std::uint64_t elapsed_ns) {
    auto& profiler_state = state();
    if (!profiler_state.enabled) {
        return;
    }

    std::scoped_lock lock(profiler_state.mutex);
    auto& phase = profiler_state.phases[std::string(phase_name)];
    ++phase.calls;
    phase.total_ns += elapsed_ns;
    phase.max_ns = std::max(phase.max_ns, elapsed_ns);
}

void set_counter(std::string_view counter_name, std::uint64_t value) {
    auto& profiler_state = state();
    if (!profiler_state.enabled) {
        return;
    }

    std::scoped_lock lock(profiler_state.mutex);
    profiler_state.counters[std::string(counter_name)] = value;
}

void write_to_configured_path() {
    auto& profiler_state = state();
    if (!profiler_state.enabled || profiler_state.output_path.empty()) {
        return;
    }

    std::unordered_map<std::string, PhaseSnapshot> phases;
    std::unordered_map<std::string, std::uint64_t> counters;
    {
        std::scoped_lock lock(profiler_state.mutex);
        phases = profiler_state.phases;
        counters = profiler_state.counters;
    }

    std::filesystem::path path(profiler_state.output_path);
    if (path.has_parent_path()) {
        std::filesystem::create_directories(path.parent_path());
    }

    std::vector<std::string> phase_names;
    phase_names.reserve(phases.size());
    for (const auto& [name, _] : phases) {
        (void)_;
        phase_names.push_back(name);
    }
    std::sort(phase_names.begin(), phase_names.end());

    std::vector<std::string> counter_names;
    counter_names.reserve(counters.size());
    for (const auto& [name, _] : counters) {
        (void)_;
        counter_names.push_back(name);
    }
    std::sort(counter_names.begin(), counter_names.end());

    std::ofstream out(path);
    if (!out.is_open()) {
        return;
    }

    out << "{\n";
    out << "  \"impl\": \"cpp\",\n";
    out << "  \"phases\": {\n";
    for (std::size_t i = 0; i < phase_names.size(); ++i) {
        const auto& name = phase_names[i];
        const auto& phase = phases[name];
        out << "    \"" << escape_json(name) << "\": {\"calls\": " << phase.calls
            << ", \"total_ns\": " << phase.total_ns
            << ", \"max_ns\": " << phase.max_ns << "}";
        out << (i + 1 == phase_names.size() ? "\n" : ",\n");
    }
    out << "  },\n";
    out << "  \"counters\": {\n";
    for (std::size_t i = 0; i < counter_names.size(); ++i) {
        const auto& name = counter_names[i];
        out << "    \"" << escape_json(name) << "\": " << counters[name];
        out << (i + 1 == counter_names.size() ? "\n" : ",\n");
    }
    out << "  }\n";
    out << "}\n";
}

ScopedPhase::ScopedPhase(std::string_view phase_name)
    : phase_name_(phase_name),
      start_(std::chrono::steady_clock::now()),
      active_(is_enabled()) {}

ScopedPhase::~ScopedPhase() {
    if (!active_) {
        return;
    }
    const auto elapsed =
        std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start_).count();
    add_duration(phase_name_, static_cast<std::uint64_t>(elapsed));
}

}
