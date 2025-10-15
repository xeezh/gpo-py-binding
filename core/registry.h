#pragma once
#include "api.h"
#include <unordered_map>
#include <vector>
#include <optional> 

namespace core {
    class Registry {
    public:
        static Registry& instance();

        void register_algo(const AlgoInfo& info);
        const AlgoInfo& get(const std::string& key) const;
        std::vector<AlgoInfo> list(std::optional<std::string> by_problem = std::nullopt) const;

    private:
        std::unordered_map<std::string, AlgoInfo> map_;
    };

#define REGISTER_ALGO(KEY, PROBLEM, SCHEMA_JSON, FN) \
    namespace {                                           \
        struct AutoReg_##KEY {                            \
            AutoReg_##KEY() {                             \
                core::Registry::instance().register_algo(  \
                    core::AlgoInfo{#KEY, PROBLEM, SCHEMA_JSON, FN}); \
            }                                             \
        } AutoRegInstance_##KEY;                          \
    }
} 
