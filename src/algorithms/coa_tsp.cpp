// cuckoo_tsp.cpp
#include "../../include/coa_tsp.h"
#include <iostream>
#include <algorithm>
#include <numeric>
#include <cmath>

using namespace coa_tsp;

    CuckooSearch::CuckooSearch(my_matrix::Matrix& distance_matrix,size_t population_size,double pa,double step_size,int max_iterations):
          distance_matrix_(distance_matrix),population_size_(population_size),pa_(pa),step_size_(step_size),max_iterations_(max_iterations),best_distance_(std::numeric_limits<double>::max()),gen_(rd_()) {
        if (distance_matrix_.GetRows() == 0 || distance_matrix_.GetRows() != distance_matrix_.GetCols()) {
            throw std::invalid_argument("Distance matrix must be square and non-empty");
        }
    }

    std::vector<size_t> CuckooSearch::solveTSP() {
    const size_t num_cities = distance_matrix_.GetRows();

    // Инициализация популяции
    std::vector<std::vector<size_t>> population(population_size_);
    std::vector<double> fitness(population_size_);

    for (size_t i = 0; i < population_size_; ++i) {
        population[i] = generateRandomSolution();
        fitness[i] = calculateDistance(population[i]);

        // Обновление лучшего решения
        if (fitness[i] < best_distance_) {
            best_distance_ = fitness[i];
            best_solution_ = population[i];
        }
    }

    for (int iter = 0; iter < max_iterations_; ++iter) {
        // Генерация нового решения с помощью Леви-полета
        for (size_t i = 0; i < population_size_; ++i) { // Levy flight for each cuckoo
            auto new_solution = levyFlight(population[i]);
            double new_fitness = calculateDistance(new_solution);

            // Выбор случайного гнезда для замены
            size_t replace_idx = std::uniform_int_distribution<size_t>(0, population_size_ - 1)(gen_);

            // Если новое решение лучше, заменяем
            if (new_fitness < fitness[replace_idx]) {
                population[replace_idx] = new_solution;
                fitness[replace_idx] = new_fitness;

                // Обновление лучшего решения
                if (new_fitness < best_distance_) {
                    best_distance_ = new_fitness;
                    best_solution_ = new_solution;
                }
            }
        }

        // Обнаружение и отбрасывание худших решений с вероятностью pa
        for (size_t i = 0; i < population_size_; ++i) {
            if (std::uniform_real_distribution<double>(0, 1)(gen_) < pa_) {
                population[i] = generateRandomSolution();
                fitness[i] = calculateDistance(population[i]);

                // Обновление лучшего решения
                if (fitness[i] < best_distance_) {
                    best_distance_ = fitness[i];
                    best_solution_ = population[i];
                }
            }
        }
       // std::cout << "Iteration " << iter << ": Best distance = " << best_distance_ << std::endl; // Всегда выводим
    }

    return best_solution_;
}


    std::vector<size_t> CuckooSearch::generateRandomSolution() {
        const size_t num_cities = distance_matrix_.GetRows();
        std::vector<size_t> solution(num_cities);

        // Заполнение последовательностью городов
        std::iota(solution.begin(), solution.end(), 0);

        // Перемешивание (кроме первого города)
        std::shuffle(solution.begin(), solution.end(), gen_); // Shuffle all cities

        return solution;
    }

double CuckooSearch::calculateDistance(const std::vector<size_t>& solution) {
        double distance = 0.0;
        const size_t n = solution.size();

        for (size_t i = 0; i < n - 1; ++i) {
            distance += distance_matrix_[solution[i]][solution[i+1]];
        }

        // Возвращение к начальному городу
        distance += distance_matrix_[solution[n-1]][solution[0]]; // Используем n-1 вместо back()

        return distance;
    }


std::vector<size_t> CuckooSearch::levyFlight(const std::vector<size_t>& current_solution) {
    std::vector<size_t> new_solution = current_solution;
    size_t n = current_solution.size();

    // Choose two distinct random cities
    size_t i = std::uniform_int_distribution<size_t>(0, n - 1)(gen_);
    size_t j = std::uniform_int_distribution<size_t>(0, n - 1)(gen_);
    while (i == j) {
        j = std::uniform_int_distribution<size_t>(0, n - 1)(gen_);
    }

    // Perform 2-opt swap
    size_t start = std::min(i, j);
    size_t end = std::max(i, j);
    std::reverse(new_solution.begin() + start + 1, new_solution.begin() + end + 1);

    return new_solution;
}

//В текущем коде этот метод не используется
    std::vector<size_t> CuckooSearch::mutateSolution(const std::vector<size_t>& solution) {
        auto mutated = solution;
        size_t n = mutated.size();

        if (n > 2) {
            size_t i = std::uniform_int_distribution<size_t>(1, n - 1)(gen_);
            size_t j = std::uniform_int_distribution<size_t>(1, n - 1)(gen_);
            std::swap(mutated[i], mutated[j]);
        }

        return mutated;
    }
