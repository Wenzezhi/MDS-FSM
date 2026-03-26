#include "mds/dfs_code.hpp"

namespace mds {

bool DfsCode::operator<(const DfsCode& other) const {
    const bool forward_first = from < to;
    const bool forward_second = other.from < other.to;

    if (!forward_first && forward_second) {
        return true;
    }
    if (forward_first && !forward_second) {
        return false;
    }

    if (forward_first && forward_second) {
        if (from != other.from) {
            return other.from < from;
        }
        if (from_label != other.from_label) {
            return from_label < other.from_label;
        }
        if (edge_label != other.edge_label) {
            return edge_label < other.edge_label;
        }
        return to_label < other.to_label;
    }

    if (to != other.to) {
        return to < other.to;
    }
    return edge_label < other.edge_label;
}

std::int32_t build_right_most_path(const std::vector<DfsCode>& dfs_codes,
                                   std::vector<std::int32_t>& right_most_path) {
    right_most_path.clear();
    std::int32_t prev_id = -1;
    for (std::size_t i = dfs_codes.size(); i-- > 0;) {
        const auto& code = dfs_codes[i];
        if (code.from < code.to && (right_most_path.empty() || prev_id == code.to)) {
            prev_id = code.from;
            right_most_path.push_back(code.to);
        }
    }
    right_most_path.push_back(0);
    return static_cast<std::int32_t>(right_most_path.size());
}

}
