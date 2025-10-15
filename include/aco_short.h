#ifndef ACO_SHORT
#define ACO_SHORT

#include "Matrix.h"
#include <vector>
#include <memory>

using namespace my_matrix;

namespace aco_short_way {
    struct AntPath {
        std::vector<std::size_t> vertices;
        double distance = 0;
    };

    /* Муравей */
    struct Ant {
        // Конструируется от стартовой вершины и конечной точки
        Ant(std::size_t start_vertex = 0, std::size_t end_vertex = 0)
            : start_location(start_vertex), current_location(start_vertex), end_location(end_vertex) {}

        // Пройденный путь и его длина
        AntPath path;

        // Пройденные вершины
        std::size_t visited_mask = 0;

        // Начальная, текущая и конечная вершины
        std::size_t start_location = 0, current_location = 0, end_location = 0;
        bool can_continue = true;

        void MakeChoice(const Matrix& g, const Matrix& p, double a, double b);

        static double getRandomChoice();
        [[nodiscard]] std::vector<std::size_t> getNeighborVertexes(const Matrix& g) const;
    };

    /* Колония муравьев */
    class AntColonyOptimization {
    public:
        AntColonyOptimization(Matrix& graph, std::size_t start, std::size_t end);
        ~AntColonyOptimization() = default;

        // Возвращает самый выгодный путь
        AntPath SolveShortestPath();

    private:
        // Определяются и подбираются опытным путем
        const double kAlpha_ = 1.0;         // Коэффициент важности феромонов (Возведение в степень)
        const double kBeta_ = 6.0;          // Коэффициент важности расстояния (Возведение в степень)
        const double kPheromone0_ = 1.0;      // Начальное количество феромонов на ребрах графа
        const double kQ_ = 5.0;             // Чем выше значение, тем больше феромонов будет откладываться
        const double kEvaporation_ = 0.0001;   // Коэффициент испарения феромонов

        Matrix* graph_;
        std::unique_ptr<Matrix> pheromone_;

        std::vector<Ant> ants_;

        // Начальная и конечная точки
        std::size_t start_location_;
        std::size_t end_location_;

        // Создает муравьев и задает начальную инициализацию каждому муравью
        void CreateAnts();

        // Глобальное обновление феромона, принимая на вход локальную составную обновления феромона
        void UpdateGlobalPheromone(const Matrix& lpu) const;
    };
}

#endif