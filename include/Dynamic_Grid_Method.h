// Dynamic_Grid_Method.h
#ifndef DYNAMIC_GRID_METHOD_H
#define DYNAMIC_GRID_METHOD_H

#include "Matrix.h"
#include <vector>
#include <limits>
#include <random>

namespace dynamic_grid {

    struct Solution {
        std::vector<size_t> path;
        double distance;

        Solution();
        void evaluate(const my_matrix::Matrix& dist_matrix);
        bool operator<(const Solution& other) const;
        bool operator==(const Solution& other) const;
    };

    struct GridParams {
        size_t initial_population = 100;
        size_t max_iterations = 1000;
        double local_expansion_factor = 0.5;
        double global_expansion_factor = 0.2;
        double min_distance = 0.1;
        bool use_3opt = false;
        bool parallel = true;
        size_t max_local_iter = 10;
        bool verbose = true;
        size_t print_interval = 10;
    };

    class DynamicGridSolver {
    public:
        DynamicGridSolver(const my_matrix::Matrix& dist_matrix, const GridParams& params);
        Solution solve();

    private:
        void validate_inputs();
        std::vector<Solution> initialize_population();
        void expand_grid(std::vector<Solution>& population);
        void mutate(Solution& sol);
        void reduce_grid(std::vector<Solution>& population);
        double solution_distance(const Solution& a, const Solution& b);
        void local_search(Solution& sol);
        bool try_3opt(Solution& sol);

        const my_matrix::Matrix& dist_matrix_;
        GridParams params_;
        std::mt19937 gen_;
        std::uniform_real_distribution<double> dist_;
    };
}
#endif
