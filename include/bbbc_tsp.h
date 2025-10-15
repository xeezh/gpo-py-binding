#ifndef BBBC_TSP
#define BBBC_TSP

#include "Matrix.h"
#include <vector>
#include <random>
#include <chrono>

using namespace my_matrix;

namespace bbbc {

    class BigBangBigCrunch {
    private:
        using Route = std::vector<int>;

        // Параметры алгоритма
        int Population_Size;        // Размер популяции
        int Max_Iterations;         // Максимальное число итераций
        int Elite_Size;             // Размер элитного пула
        double Mutation_Rate;       // Вероятность мутации
        size_t Max_Neighbors;  // Максимальное число соседей для локального поиска


        Matrix& graph_;             // Матрица расстояний между городами

        std::vector<Route> population;     // Текущая популяция маршрутов
        Route bestRoute;                   // Лучший найденный маршрут
        double bestDistance;              // Его длина

        std::vector<Route> elitePool;      // Элитный пул (хорошие решения)
        std::vector<double> eliteDistances;// Длины маршрутов в элитном пуле

        std::default_random_engine rng;    // Генератор случайных чисел

        // Генерация начальной популяции случайных маршрутов
        void initializePopulation();

        // Расчёт длины маршрута (сумма расстояний между городами по кругу)
        double Tour_Distance_Calc(const Route& tour);

        // Генерация соседних маршрутов (например, перестановками)
        std::vector<Route> createNeighbors(const Route& route ) const;

        // Локальный поиск вокруг текущего маршрута (например, 2-opt)
        double localSearch(Route& route,size_t max_neighbors);
        bool isValidTour(const std::vector<int>& tour);
        // Обновление элитного пула новыми маршрутами
        void updateElitePool(const Route& newRoute, double newDistance);

        // Расчёт "центра масс" из элитных решений
        Route calculateCenterOfMass();

        // Удаление худших решений из популяции
        void removeWorstSolutions();

    public:

        BigBangBigCrunch(Matrix& graph,
                         int populationSize,
                         int maxIterations,
                         int eliteSize,
                         double mutationRate,
                         unsigned int seed = std::random_device{}(),
                         size_t Max_Neighbors = 50 );

        ~BigBangBigCrunch();

        // Получить длину лучшего найденного маршрута
        double getBestDistance() const;

        // Получить сам лучший маршрут
        const Route& getBestRoute() const;

        // Запуск алгоритма решения задачи коммивояжёра
        void solve();

        // Напечатать лучший маршрут и его длину
        void printBestRoute() const;
    };

}

#endif
