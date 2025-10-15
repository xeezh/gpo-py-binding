#include "../../include/Greedy_Adaptive_Method.h"
#include <numeric>
#include <algorithm>
#include <iostream>
#include <climits>

namespace greedy_adaptive {

void Solution::evaluate(const my_matrix::Matrix& dist_matrix) {
    distance = 0.0;
    for(size_t i = 0; i < path.size() - 1; ++i) {
        distance += dist_matrix[path[i]][path[i+1]];
    }
    distance += dist_matrix[path.back()][path[0]];
}

GreedyAdaptiveSolver::GreedyAdaptiveSolver(const my_matrix::Matrix& dist_matrix, size_t iterations)
    : dist_matrix_(dist_matrix),
      max_iterations_(iterations),
      rcl_min_size_(3),
      rcl_max_size_(15),
      elite_pool_size_(5)
{
    std::random_device rd;
    gen_.seed(rd());
}

std::vector<size_t> GreedyAdaptiveSolver::constructInitialSolution() {
    std::vector<size_t> candidates(dist_matrix_.GetRows());
    std::iota(candidates.begin(), candidates.end(), 0);

    std::vector<size_t> path;
    path.reserve(dist_matrix_.GetRows());

    // Адаптивный расчет RCL
    auto calculateRCLSize = [this](size_t candidates_left) -> size_t {
        return std::clamp(candidates_left / 4, rcl_min_size_, rcl_max_size_);
    };

    // Выбор стартового города
    std::uniform_int_distribution<size_t> dist(0, candidates.size()-1);
    size_t current = dist(gen_);
    path.push_back(candidates[current]);
    candidates.erase(candidates.begin() + current);

    while(!candidates.empty()) {
        std::vector<std::pair<size_t, double>> evaluated;
        for(size_t i = 0; i < candidates.size(); ++i) {
            evaluated.emplace_back(i, dist_matrix_[path.back()][candidates[i]]);
        }

        // Сортировка и выбор из RCL
        std::sort(evaluated.begin(), evaluated.end(),
            [](const auto& a, const auto& b) { return a.second < b.second; });

        size_t RCL_size = calculateRCLSize(candidates.size());
        RCL_size = std::min(RCL_size, evaluated.size());

        std::uniform_int_distribution<size_t> rnd_dist(0, RCL_size-1);
        size_t selected_idx = evaluated[rnd_dist(gen_)].first;

        path.push_back(candidates[selected_idx]);
        candidates.erase(candidates.begin() + selected_idx);
    }
    return path;
}

std::vector<size_t> GreedyAdaptiveSolver::twoOptSwap(
    const std::vector<size_t>& path,
    size_t i,
    size_t k) const
{
    std::vector<size_t> new_path(path.begin(), path.begin() + i);
    std::reverse_copy(path.begin() + i, path.begin() + k + 1, back_inserter(new_path));
    new_path.insert(new_path.end(), path.begin() + k + 1, path.end());
    return new_path;
}

std::vector<size_t> GreedyAdaptiveSolver::localSearch(const std::vector<size_t>& path) {
    std::vector<size_t> best_path = path;
    double best_dist = calculateDistance(path);
    bool improvement = true;

    // Полный 2-opt поиск
    while(improvement) {
        improvement = false;

        for(size_t i = 1; i < path.size()-1; ++i) {
            for(size_t k = i+1; k < path.size(); ++k) {
                auto new_path = twoOptSwap(best_path, i, k);
                double new_dist = calculateDistance(new_path);

                if(new_dist < best_dist) {
                    best_path = new_path;
                    best_dist = new_dist;
                    improvement = true;
                }
            }
        }
    }
    return best_path;
}

double GreedyAdaptiveSolver::calculateDistance(const std::vector<size_t>& path) const {
    double total = 0.0;
    for(size_t i = 0; i < path.size()-1; ++i) {
        total += dist_matrix_[path[i]][path[i+1]];
    }
    return total + dist_matrix_[path.back()][path[0]];
}

Solution GreedyAdaptiveSolver::solve() {
    Solution best_solution;
    best_solution.distance = std::numeric_limits<double>::max();

    std::vector<Solution> elite_pool;
    elite_pool.reserve(elite_pool_size_);

    for(size_t iter = 0; iter < max_iterations_; ++iter) {
        auto current_path = constructInitialSolution();
        current_path = localSearch(current_path);

        Solution current_solution;
        current_solution.path = current_path;
        current_solution.evaluate(dist_matrix_);

        // Обновление элитного пула
        if(elite_pool.size() < elite_pool_size_) {
            elite_pool.push_back(current_solution);
        } else {
            auto worst = std::max_element(elite_pool.begin(), elite_pool.end(),
                [](const auto& a, const auto& b) { return a.distance < b.distance; });

            if(current_solution.distance < worst->distance) {
                *worst = current_solution;
            }
        }

        // Обновление глобального лучшего
        if(current_solution.distance < best_solution.distance) {
            best_solution = current_solution;
            //std::cout << "Iter " << iter << " | Best: " << best_solution.distance << "\n";
        }
    }

    // Финальная оптимизация лучшего решения
    best_solution.path = localSearch(best_solution.path);
    best_solution.evaluate(dist_matrix_);

    return best_solution;
}

} // namespace greedy_adaptive
