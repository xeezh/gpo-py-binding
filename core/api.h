#pragma once
#include <vector>
#include <string>
#include <functional>
#include <optional>
#include <stdexcept>
#include <nlohmann/json.hpp>
#include "../include/Matrix.h"

namespace core {

    struct Result {
        std::vector<int> path;      
        double cost = 0.0;          
        std::string algorithm;      
        std::string problem;        
        nlohmann::json meta;        
    };

    using Params = nlohmann::json;

    using AlgoFn = std::function<Result(const my_matrix::Matrix&, const Params&)>;

    struct AlgoInfo {
        std::string key;           
        std::string problem;     
        nlohmann::json schema;  
        AlgoFn run;
    };

    struct UavError : std::runtime_error {
        using std::runtime_error::runtime_error;
    };

}