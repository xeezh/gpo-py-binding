#ifndef ACO_TSP
#define ACO_TSP

#include "Matrix.h"
#include "vector"
#include <memory>

using namespace my_matrix;

namespace aco_tsp {
    /* Муравей */
    struct Ant {
        // Конструктор по умолчанию (если нужно)
        Ant() : start_location(0), current_location(0), can_continue(true), visited_mask(0) {}
        // Конструируется от стартовой вершины
        explicit  Ant(std::size_t start) : start_location(start), current_location(start), can_continue(true), visited_mask(0) {
            visited_mask |= (1ULL << start_location); // Отмечаем начальную вершину как посещенную
        }

        // Пройденные вершины
        std::size_t visited_mask = 0;

        // начальная и текущая вершины
        std::size_t start_location = 0, current_location = 0;

        // Флаг для возможности продолжения пути
        bool can_continue = true;

        // Выбор следующей вершины
        void MakeChoice(const Matrix& g, const Matrix& p, double a, double b);

        // Генерация случайного числа
        static double getRandomChoice();

        // Получение соседей текущей вершины
        [[nodiscard]] std::vector<std::size_t> getNeighborVertexes(const Matrix& g) const;
    };

    /* Путь муравья */
    struct AntPath {
        std::vector<std::size_t> vertices;
        double distance;

        // Конструктор по умолчанию
        AntPath() : distance(0.0) {}

        // Конструктор с параметрами
        AntPath(const std::vector<std::size_t>& vertices, double distance)
                : vertices(vertices), distance(distance) {}
    };

    /* Колония муравьев */
    class AntColonyOptimization {
    public:
        explicit AntColonyOptimization(Matrix& graph, std::size_t evaporation, std::size_t initial_pheromone, std::size_t q);
        ~AntColonyOptimization() = default;

        // Сброс феромонов
        void ResetPheromones();

        // Решение задачи коммивояжера
        AntPath SolveSalesmansProblem(std::size_t max_iterations,
                             std::size_t max_iterations_without_improvement,
                             std::size_t elite_ants_count);

    private:
        double kAlpha_;             // Коэффициент важности феромонов
        double kBeta_;              // Коэффициент важности расстояния
        double kEvaporation_;       // Коэффициент испарения феромонов
        const double kPheromone0_;  // Начальное количество феромонов на ребрах
        const double kQ_;           // Параметр для расчета феромонов

        Matrix* graph_;
        std::unique_ptr<Matrix> pheromone_;

        std::vector<Ant> ants_;  // Муравьи в колонии

        // Создание муравьев
        void CreateAnts();

        // Обновление глобальных феромонов
        void UpdateGlobalPheromone(const Matrix& lpu, const std::vector<AntPath>& elite_paths) const;

        // Обновление локальных феромонов
        void UpdateLocalPheromone(std::size_t from, std::size_t to);

        // Регулировка параметров по мере выполнения
        void AdjustParameters(std::size_t iteration, std::size_t max_iterations);
    };
}

#endif