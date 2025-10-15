// Dynamic_Grid_Method.cpp
#include "../../include/Dynamic_Grid_Method.h"
#include <algorithm>
#include <numeric>
#include <random>
#include <stdexcept>
#include <unordered_set>
#include <iostream>
#include <chrono>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace dynamic_grid {

    Solution::Solution() : distance(std::numeric_limits<double>::max()) {}

    void Solution::evaluate(const my_matrix::Matrix& dist_matrix) {
        distance = 0.0;
        if (path.empty()) return;

        for (size_t i = 0; i < path.size(); ++i) {
            size_t next = (i + 1) % path.size();
            distance += dist_matrix[path[i]][path[next]];
        }
    }

    bool Solution::operator<(const Solution& other) const {
        return distance < other.distance;
    }
    bool Solution::operator==(const Solution& other) const
    {
        return path == other.path;
    }

    DynamicGridSolver::DynamicGridSolver(const my_matrix::Matrix& dist_matrix, const GridParams& params)
        : dist_matrix_(dist_matrix), params_(params), gen_(std::random_device{}()), dist_(0.0, 1.0) {
        validate_inputs();
    }

    void DynamicGridSolver::validate_inputs() {
        if (dist_matrix_.GetRows() == 0 || dist_matrix_.GetCols() == 0) {
            throw std::invalid_argument("Distance matrix cannot be empty");
        }
        if (params_.initial_population == 0) {
            throw std::invalid_argument("Initial population size cannot be zero");
        }
    }

Solution DynamicGridSolver::solve() {
    std::vector<Solution> population = initialize_population();
    if (population.empty()) return Solution();

    Solution global_best = *std::min_element(population.begin(), population.end());


    auto start_time = std::chrono::high_resolution_clock::now();
    size_t last_printed_iter = 0;

    for (size_t iter = 0; iter < params_.max_iterations; ++iter) {
        expand_grid(population);

        if (params_.parallel) {
            #pragma omp parallel for
            for (size_t i = 0; i < population.size(); ++i) {
                local_search(population[i]);
            }
        } else {
            for (auto& sol : population) {
                local_search(sol);
            }
        }

        reduce_grid(population);

        if (!population.empty()) {
            auto current_best = *std::min_element(population.begin(), population.end());
            if (current_best < global_best) {
                global_best = current_best;

                // Выводим только при улучшении и с заданным интервалом
                if (params_.verbose && (iter - last_printed_iter >= params_.print_interval || iter == 0)) {
                    auto current_time = std::chrono::high_resolution_clock::now();
                    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
                        current_time - start_time).count() / 1000.0;


                    last_printed_iter = iter;
                }
            }
        }
    }

    if (params_.verbose) {
        auto end_time = std::chrono::high_resolution_clock::now();
        auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>(
            end_time - start_time).count() / 1000.0;



    }

    return global_best;
}

    std::vector<Solution> DynamicGridSolver::initialize_population() {
        std::vector<Solution> population;
        size_t n = dist_matrix_.GetRows();
        if (n == 0) return population;

        population.reserve(params_.initial_population);
        std::vector<size_t> base_path(n);
        std::iota(base_path.begin(), base_path.end(), 0);

        std::random_device rd;
        std::mt19937 gen(rd());

        for (size_t i = 0; i < params_.initial_population; ++i) {
            Solution sol;
            sol.path = base_path;
            std::shuffle(sol.path.begin(), sol.path.end(), gen);
            sol.evaluate(dist_matrix_);
            population.push_back(sol);
        }

        return population;
    }

    void DynamicGridSolver::expand_grid(std::vector<Solution>& population) {
        const size_t current_size = population.size();
        if (current_size == 0) return;

        std::vector<Solution> new_solutions;
        size_t local_new = static_cast<size_t>(current_size * params_.local_expansion_factor);
        size_t global_new = static_cast<size_t>(current_size * params_.global_expansion_factor);
        new_solutions.reserve(local_new + global_new);

        // Local expansion with mutation
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<size_t> idx_dist(0, current_size - 1);
        std::uniform_real_distribution<double> real_dist(0.0, 1.0);

        for (size_t i = 0; i < local_new; ++i) {
            Solution new_sol = population[idx_dist(gen)];
            mutate(new_sol);
            new_sol.evaluate(dist_matrix_);
            new_solutions.push_back(new_sol);
        }

        // Global expansion with new random solutions
        size_t n = dist_matrix_.GetRows();
        std::vector<size_t> base_path(n);
        std::iota(base_path.begin(), base_path.end(), 0);

        for (size_t i = 0; i < global_new; ++i) {
            Solution new_sol;
            new_sol.path = base_path;
            std::shuffle(new_sol.path.begin(), new_sol.path.end(), gen);
            new_sol.evaluate(dist_matrix_);
            new_solutions.push_back(new_sol);
        }

        population.insert(population.end(), new_solutions.begin(), new_solutions.end());
    }

    void DynamicGridSolver::mutate(Solution& sol) {
        if (sol.path.size() < 2) return;

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<size_t> idx_dist(0, sol.path.size() - 1);
        std::uniform_real_distribution<double> real_dist(0.0, 1.0);

        size_t i = idx_dist(gen);
        size_t j = idx_dist(gen);

        if (i > j) std::swap(i, j);

        if (real_dist(gen) < 0.5) {
            std::reverse(sol.path.begin() + i, sol.path.begin() + j + 1);
        }
        else {
            std::swap(sol.path[i], sol.path[j]);
        }
    }

    void DynamicGridSolver::reduce_grid(std::vector<Solution>& population) {
        if (population.size() <= 1) return;

        // Remove duplicates
        std::unordered_set<std::string> seen;
        auto new_end = std::remove_if(population.begin(), population.end(),
            [&seen](const Solution& s) {
                std::string key;
                key.reserve(s.path.size() * sizeof(size_t));
                for (auto node : s.path) {
                    key.append(reinterpret_cast<const char*>(&node), sizeof(node));
                }
                return !seen.insert(key).second;
            });
        population.erase(new_end, population.end());

        // Sort and keep best solutions
        std::sort(population.begin(), population.end());
        if (population.size() > params_.initial_population) {
            population.resize(params_.initial_population);
        }
    }

    void DynamicGridSolver::local_search(Solution& sol) {
        const size_t n = sol.path.size();
        if (n < 4) return;

        bool improved = true;
        size_t iteration_count = 0;

        while (improved && iteration_count < params_.max_local_iter) {
            improved = false;
            iteration_count++;

            // 2-opt optimization
            for (size_t i = 0; i < n - 1 && !improved; ++i) {
                size_t a = sol.path[i];
                size_t b = sol.path[i + 1];
                double ab = dist_matrix_[a][b];

                for (size_t j = i + 2; j < n && !improved; ++j) {
                    size_t c = sol.path[j];
                    size_t d = sol.path[(j + 1) % n];
                    double cd = dist_matrix_[c][d];
                    double current = ab + cd;

                    double ac = dist_matrix_[a][c];
                    double bd = dist_matrix_[b][d];
                    double new_dist = ac + bd;

                    if (new_dist < current - 1e-9) {
                        std::reverse(sol.path.begin() + i + 1, sol.path.begin() + j + 1);
                        sol.distance += (new_dist - current);
                        improved = true;
                    }
                }
            }

            // Optional 3-opt
            if (!improved && params_.use_3opt && n >= 6) {
                improved = try_3opt(sol);
            }
        }
    }

    bool DynamicGridSolver::try_3opt(Solution& sol) {
        const size_t n = sol.path.size();
        bool improved = false;

        for (size_t i = 1; i < n - 4 && !improved; ++i) {
            for (size_t j = i + 2; j < n - 2 && !improved; ++j) {
                for (size_t k = j + 2; k < n && !improved; ++k) {
                    size_t a = sol.path[i - 1];
                    size_t b = sol.path[i];
                    size_t c = sol.path[j - 1];
                    size_t d = sol.path[j];
                    size_t e = sol.path[k - 1];
                    size_t f = sol.path[k % n];

                    double current = dist_matrix_[a][b] + dist_matrix_[c][d] + dist_matrix_[e][f];
                    double new_dist = dist_matrix_[a][d] + dist_matrix_[b][f] + dist_matrix_[c][e];

                    if (new_dist < current - 1e-9) {
                        std::reverse(sol.path.begin() + i, sol.path.begin() + j);
                        std::reverse(sol.path.begin() + j, sol.path.begin() + k);
                        sol.distance += (new_dist - current);
                        improved = true;
                    }
                }
            }
        }

        return improved;
    }
}
