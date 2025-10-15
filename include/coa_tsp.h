// cuckoo_tsp.h
#ifndef COA_TSP_H
#define COA_TSP_H

#include "Matrix.h"
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>

using namespace my_matrix;;

namespace coa_tsp {

    class CuckooSearch {
    public:
        explicit CuckooSearch(my_matrix::Matrix& distance_matrix, // матрицу расстояний между городами (квадратную)
                              size_t population_size = 25,  //размер популяции решений (количество "гнёзд").
                              double pa = 0.25,  //вероятность отбрасывания худших решений (параметр алгоритма).
                              double step_size = 1.0,  //размер шага для Леви-полёта (здесь не используется явно).
                              int max_iterations = 1000);// максимальное число итераций алгоритма. сколько раз алгоритм будет пытаться улучшить решения.

        // Основная функция решения
        std::vector<size_t> solveTSP();

        double getBestDistance() const { return best_distance_; }

    private:
        std::vector<size_t> generateRandomSolution();

        double calculateDistance(const std::vector<size_t>& solution);

        // Леви-полет для генерации нового решения
        std::vector<size_t> levyFlight(const std::vector<size_t>& current_solution);

        // Мутация решения (обмен двух городов)
        std::vector<size_t> mutateSolution(const std::vector<size_t>& solution);

        my_matrix::Matrix& distance_matrix_;
        size_t population_size_;
        double pa_; // Вероятность обнаружения
        double step_size_;
        int max_iterations_;
        double best_distance_;
        std::vector<size_t> best_solution_;

        std::random_device rd_;
        std::mt19937 gen_;

        //void twoOpt(std::vector<size_t>& route);
    };

}

#endif