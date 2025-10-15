#include "registry.h"

using namespace core;

Registry& Registry::instance() {
    static Registry R;
    return R;
}

void Registry::register_algo(const AlgoInfo& info) {
    map_[info.key] = info;
}

const AlgoInfo& Registry::get(const std::string& key) const {
    auto it = map_.find(key);
    if (it == map_.end()) throw UavError("Unknown algorithm: " + key);
    return it->second;
}

std::vector<AlgoInfo> Registry::list(std::optional<std::string> by_problem) const {
    std::vector<AlgoInfo> out;
    out.reserve(map_.size());
    for (auto& [k, v] : map_) {
        if (!by_problem || v.problem == *by_problem) out.push_back(v);
    }
    return out;
}