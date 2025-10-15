#ifndef GREEDY_ADAPTIVE_METHOD_H
#define GREEDY_ADAPTIVE_METHOD_H

#include "Matrix.h"
#include <vector>
#include <random>
#include <algorithm>
#include <limits>  // Для numeric_limits

namespace greedy_adaptive {

    struct Solution {
        std::vector<size_t> path;
        double distance;

        Solution() : distance(std::numeric_limits<double>::max()) {}
        void evaluate(const my_matrix::Matrix& dist_matrix);
    };

    class GreedyAdaptiveSolver {
    public:
        GreedyAdaptiveSolver(const my_matrix::Matrix& dist_matrix, size_t iterations);
        Solution solve();

    private:
        // Основные методы
        std::vector<size_t> constructInitialSolution();
        std::vector<size_t> localSearch(const std::vector<size_t>& path);
        std::vector<size_t> twoOptSwap(const std::vector<size_t>& path, size_t i, size_t k) const;
        double calculateDistance(const std::vector<size_t>& path) const;

        // Вспомогательные параметры
        size_t rcl_min_size_ = 3;    // Минимальный размер RCL
        size_t rcl_max_size_ = 15;   // Максимальный размер RCL
        size_t elite_pool_size_ = 5; // Размер элитного пула

        // Основные параметры
        const my_matrix::Matrix& dist_matrix_;
        size_t max_iterations_;
        std::mt19937 gen_;
    };

} // namespace greedy_adaptive

#endif // GREEDY_ADAPTIVE_METHOD_H
