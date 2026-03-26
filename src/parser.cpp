#include "mds/parser.hpp"

#include <algorithm>
#include <charconv>
#include <fstream>
#include <stdexcept>
#include <string_view>
#include <tuple>
#include <vector>

namespace mds {

namespace {

std::int32_t parse_int32(std::string_view text, const char* error_message) {
    std::int32_t value = 0;
    const auto* begin = text.data();
    const auto* end = begin + text.size();
    const auto [ptr, ec] = std::from_chars(begin, end, value);
    if (ec != std::errc{} || ptr != end) {
        throw std::runtime_error(error_message);
    }
    return value;
}

std::int32_t parse_next_int32(const char*& cursor, const char* end, const char* error_message) {
    while (cursor < end && *cursor == ' ') {
        ++cursor;
    }
    std::int32_t value = 0;
    const auto [ptr, ec] = std::from_chars(cursor, end, value);
    if (ec != std::errc{}) {
        throw std::runtime_error(error_message);
    }
    cursor = ptr;
    return value;
}

}

ParseResult read_gfu_format(const std::string& file_path) {
    std::ifstream input(file_path, std::ios::binary | std::ios::ate);
    if (!input.is_open()) {
        throw std::runtime_error("Failed to open file: " + file_path);
    }

    const auto file_size = input.tellg();
    if (file_size < 0) {
        throw std::runtime_error("Failed to read file size: " + file_path);
    }
    std::string file_content(static_cast<std::size_t>(file_size), '\0');
    input.seekg(0, std::ios::beg);
    input.read(file_content.data(), file_size);
    if (!input && !input.eof()) {
        throw std::runtime_error("Failed to read file: " + file_path);
    }

    ParseResult result;
    std::vector<std::tuple<std::int32_t, std::int32_t, std::int32_t>> edges;
    std::string_view remaining(file_content);

    auto next_line = [&](const char* error_message) -> std::string_view {
        if (remaining.empty()) {
            throw std::runtime_error(error_message);
        }

        const auto line_end = remaining.find_first_of("\r\n");
        if (line_end == std::string_view::npos) {
            const auto line = remaining;
            remaining = {};
            return line;
        }

        const auto line = remaining.substr(0, line_end);
        auto next = line_end;
        while (next < remaining.size() && (remaining[next] == '\r' || remaining[next] == '\n')) {
            ++next;
        }
        remaining.remove_prefix(next);
        return line;
    };

    (void)next_line("Expected header line");

    const auto n_vertex = parse_int32(next_line("Expected number of vertices"), "Failed to parse vertex count");
    result.graph.init_vertex(n_vertex);
    result.label_map.reserve(static_cast<std::size_t>(n_vertex));
    result.label_map_inverse.reserve(static_cast<std::size_t>(n_vertex));
    result.graph.label_frequency.reserve(static_cast<std::size_t>(n_vertex));

    for (std::int32_t i = 0; i < n_vertex; ++i) {
        const auto label_line = next_line("Expected vertex label");
        const auto [it, inserted] =
            result.label_map.try_emplace(std::string(label_line), static_cast<std::int32_t>(result.label_map.size()));
        const auto label_id = it->second;
        if (inserted) {
            result.label_map_inverse.emplace(label_id, it->first);
        }
        result.graph.label[static_cast<std::size_t>(i)] = label_id;
        ++result.graph.label_frequency[label_id];
    }

    const auto n_edge = parse_int32(next_line("Expected number of edges"), "Failed to parse edge count");
    result.graph.init_edge(n_edge);
    edges.reserve(static_cast<std::size_t>(n_edge));
    result.edge_label_map.reserve(static_cast<std::size_t>(n_edge));
    result.edge_label_map_inverse.reserve(static_cast<std::size_t>(n_edge));
    result.edge_list.reserve(static_cast<std::size_t>(n_edge));

    for (std::int32_t i = 0; i < n_edge; ++i) {
        const auto edge_line = next_line("Expected edge");
        const auto* cursor = edge_line.data();
        const auto* end = cursor + edge_line.size();
        const auto fst = parse_next_int32(cursor, end, "Failed to parse source vertex");
        const auto snd = parse_next_int32(cursor, end, "Failed to parse target vertex");
        const auto raw_label = parse_next_int32(cursor, end, "Failed to parse edge label");

        const auto [it, inserted] =
            result.edge_label_map.try_emplace(raw_label, static_cast<std::int32_t>(result.edge_label_map.size()));
        const auto edge_label_id = it->second;
        if (inserted) {
            result.edge_label_map_inverse.emplace(edge_label_id, raw_label);
        }

        edges.emplace_back(fst, snd, edge_label_id);
        const auto a = std::min(result.graph.label[static_cast<std::size_t>(fst)],
                                result.graph.label[static_cast<std::size_t>(snd)]);
        const auto b = std::max(result.graph.label[static_cast<std::size_t>(fst)],
                                result.graph.label[static_cast<std::size_t>(snd)]);
        const auto key = std::make_tuple(a, b, edge_label_id);
        if (result.edge_list.find(key) == result.edge_list.end()) {
            result.edge_list.emplace(key, static_cast<std::int32_t>(result.edge_list.size()));
        }

        result.graph.max_edge_label = std::max(result.graph.max_edge_label, edge_label_id + 1);
        ++result.graph.degree[static_cast<std::size_t>(fst)];
        ++result.graph.degree[static_cast<std::size_t>(snd)];
    }

    std::int32_t pos = 0;
    for (std::size_t i = 0; i < static_cast<std::size_t>(result.graph.n_vertex); ++i) {
        result.graph.nbr_offset[i] = pos;
        pos += result.graph.degree[i];
        result.graph.core[i] = result.graph.degree[i];
        result.graph.max_degree = std::max(result.graph.max_degree, result.graph.degree[i]);
    }
    result.graph.nbr_offset[static_cast<std::size_t>(result.graph.n_vertex)] = pos;

    std::fill(result.graph.degree.begin(), result.graph.degree.end(), 0);

    for (const auto& [fst, snd, label_id] : edges) {
        const auto fst_u = static_cast<std::size_t>(fst);
        const auto snd_u = static_cast<std::size_t>(snd);
        const auto pos_fst = static_cast<std::size_t>(result.graph.nbr_offset[fst_u] + result.graph.degree[fst_u]);
        result.graph.nbr[pos_fst] = {snd, label_id};
        ++result.graph.degree[fst_u];

        const auto pos_snd = static_cast<std::size_t>(result.graph.nbr_offset[snd_u] + result.graph.degree[snd_u]);
        result.graph.nbr[pos_snd] = {fst, label_id};
        ++result.graph.degree[snd_u];
    }

    result.graph.n_label = static_cast<std::int32_t>(result.label_map.size());
    return result;
}

}
