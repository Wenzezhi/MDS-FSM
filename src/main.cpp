#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#ifdef _DEBUG
#include <crtdbg.h>
#endif
#include <string>
#include <vector>

#include "mds/mining/core.hpp"
#include "mds/profiler.hpp"

namespace {

[[noreturn]] void print_usage() {
    std::cerr << "Usage: mds-fsm -d <dataset> -k <k> [-o <output_dir>]\n";
    std::cerr << "  -d <dataset>    Dataset file path (e.g., dataset/yeast.gfu)\n";
    std::cerr << "  -k <k>          Top-k value\n";
    std::cerr << "  -o <output_dir> Output directory (default: output/)\n";
    std::exit(1);
}

}

int main(int argc, char** argv) {
#ifdef _DEBUG
    _CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_FILE | _CRTDBG_MODE_DEBUG);
    _CrtSetReportFile(_CRT_WARN, _CRTDBG_FILE_STDERR);
    _CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_FILE | _CRTDBG_MODE_DEBUG);
    _CrtSetReportFile(_CRT_ERROR, _CRTDBG_FILE_STDERR);
    _CrtSetReportMode(_CRT_ASSERT, _CRTDBG_MODE_FILE | _CRTDBG_MODE_DEBUG);
    _CrtSetReportFile(_CRT_ASSERT, _CRTDBG_FILE_STDERR);
#endif

    std::string dataset;
    std::size_t k = 0;
    std::string output_dir = "output";

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "-d") {
            ++i;
            if (i < argc) {
                dataset = argv[i];
            }
        } else if (arg == "-k") {
            ++i;
            if (i < argc) {
                try {
                    k = static_cast<std::size_t>(std::stoul(argv[i]));
                } catch (...) {
                    k = 0;
                }
            }
        } else if (arg == "-o") {
            ++i;
            if (i < argc) {
                output_dir = argv[i];
            }
        }
    }

    if (dataset.empty() || k == 0) {
        print_usage();
    }

    const auto dataset_path =
        (dataset.find('/') != std::string::npos || dataset.find('\\') != std::string::npos ||
         dataset.ends_with(".gfu"))
            ? (dataset.ends_with(".gfu") ? dataset : dataset + ".gfu")
            : "dataset/" + dataset + ".gfu";

    const auto dataset_name =
        std::filesystem::path(dataset_path).stem().string().empty()
            ? std::string("output")
            : std::filesystem::path(dataset_path).stem().string();
    const auto output_path = output_dir + "/" + dataset_name + "_" + std::to_string(k) + ".txt";

    const auto start_time = std::chrono::steady_clock::now();

    mds::MiningContext ctx;
    {
        mds::profile::ScopedPhase phase("load_graph_total");
        ctx.load_graph(dataset_path);
    }

    std::cout << "Dataset: " << dataset_path << "\n";
    std::cout << "Vertices: " << ctx.data_graph.n_vertex << ", Edges: " << ctx.data_graph.n_edge << "\n";

    {
        mds::profile::ScopedPhase phase("process_data_graph_total");
        ctx.process_data_graph();
    }

    std::cout << "Unique vertex labels: " << ctx.n_unique_label
              << ", Unique edge labels: " << ctx.n_unique_edge_label << "\n";

    ctx.frequent_mining(k, k);

    const auto total_time =
        std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - start_time)
            .count();
    const auto total_ns =
        static_cast<std::uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(
                                        std::chrono::steady_clock::now() - start_time)
                                        .count());

    std::cout << "\nTotal runtime: " << std::fixed << std::setprecision(3) << total_time << "s\n";
    mds::profile::set_counter("total_runtime_ns", total_ns);
    mds::profile::set_counter("patterns_verified", ctx.patterns_verified);
    mds::profile::set_counter("candidates_generated", ctx.candidates_generated);
    mds::profile::set_counter("candidates_added_to_queue", ctx.candidates_added_to_queue);
    mds::profile::set_counter("answer_set_size", static_cast<std::uint64_t>(ctx.answer_set.size()));
    mds::profile::set_counter("minimality_cache_hits", ctx.minimality_cache.hits);
    mds::profile::set_counter("minimality_cache_misses", ctx.minimality_cache.misses);

    ctx.print_answer();

    std::filesystem::create_directories(output_dir);
    std::ofstream output(output_path);
    if (output.is_open()) {
        output << "Dataset: " << dataset_path << "\n";
        output << "K: " << k << "\n";
        output << "Total Time: " << std::fixed << std::setprecision(3) << total_time << "s\n\n";
        output << "Results (" << ctx.answer_set.size() << " patterns):\n";
        if (!ctx.answer_set.empty()) {
            output << "Minimum MDS: " << ctx.answer_set.front().mds << "\n";
        }
        output << "\n";

        for (std::size_t i = 0; i < ctx.answer_set.size(); ++i) {
            output << "Pattern " << (i + 1) << ": MDS = " << ctx.answer_set[i].mds << "\n";
            for (const auto& code : ctx.answer_set[i].graph.dfs_code) {
                const auto fl_it = ctx.label_map_inverse.find(code.from_label);
                const auto tl_it = ctx.label_map_inverse.find(code.to_label);
                const auto el_it = ctx.edge_label_map_inverse.find(code.edge_label);
                const auto& fl = fl_it != ctx.label_map_inverse.end() ? fl_it->second : std::string{"?"};
                const auto& tl = tl_it != ctx.label_map_inverse.end() ? tl_it->second : std::string{"?"};
                const auto el = el_it != ctx.edge_label_map_inverse.end() ? el_it->second : -1;
                output << "  ( " << code.from << ", " << code.to << ", " << fl << ", " << el << ", " << tl
                       << " )\n";
            }
        }
        std::cout << "\nResults written to: " << output_path << "\n";
    }

    mds::profile::write_to_configured_path();

    return 0;
}
