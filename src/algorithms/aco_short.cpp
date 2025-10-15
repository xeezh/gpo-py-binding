#include <iostream>
#include <chrono>
#include <atomic>
#include <random>
#include <algorithm>
#include <cassert>
#include <stdexcept>
#include "../../include/aco_short.h"


/* Конструктор колонии */
aco_short_way::AntColonyOptimization::AntColonyOptimization(Matrix& graph, std::size_t start, std::size_t end)
    : kQ_(2 * graph.GetRows()), graph_(&graph), pheromone_(nullptr), start_location_(start), end_location_(end) {
    const std::size_t kVertexesCount = graph_->GetRows();

    std::vector<double> array(kVertexesCount * kVertexesCount, kPheromone0_);
    pheromone_ = std::make_unique<Matrix>(kVertexesCount, kVertexesCount, array.data());
    CreateAnts(); // создание муравьев
}

/* Инициализация муравьев */
void aco_short_way::AntColonyOptimization::CreateAnts() {
    const auto kAntsCount = graph_->GetRows();
    ants_.resize(kAntsCount);

    for (std::size_t i = 0; i != kAntsCount; ++i) {
        ants_[i] = Ant(start_location_, end_location_); // Передаем начальную и конечную точки
    }
}

/* Метод глобального обновления феромона */
void aco_short_way::AntColonyOptimization::UpdateGlobalPheromone(const Matrix& lpu) const {
    for (std::size_t from = 0; from < lpu.GetRows(); ++from) {
        for (std::size_t to = 0; to < lpu.GetRows(); ++to) {
            // Проверка индексов
            if (from >= pheromone_->GetRows() || to >= pheromone_->GetRows())
                throw std::out_of_range("Index out of bounds in UpdateGlobalPheromone");

            (*pheromone_)[from][to] = (1 - kEvaporation_) * (*pheromone_)[from][to] + lpu[from][to];
            if ((*pheromone_)[from][to] < 0.01 && from != to) {
                (*pheromone_)[from][to] = 0.01;
            }
        }
    }
}

/* Получение случайного числа в пределах (0:1) */
double aco_short_way::Ant::getRandomChoice() {
    std::uniform_real_distribution<> dist(0.0, 1.0);
    std::default_random_engine re(std::chrono::system_clock::now().time_since_epoch().count());
    return dist(re);
}

std::vector<std::size_t> aco_short_way::Ant::getNeighborVertexes(const Matrix& g) const {
    std::vector<std::size_t> vertexes;
    for (std::size_t to = 0; to != g.GetRows(); ++to) {
        if (current_location >= g.GetRows())
            throw std::out_of_range("current_location is out of bounds");

        bool edge_is_exist = g[current_location][to] != 0;
        bool vertex_is_unvisited = !(visited_mask & (1ULL << to));
        if (edge_is_exist && vertex_is_unvisited)
            vertexes.push_back(to);
    }
    return vertexes;
}

void aco_short_way::Ant::MakeChoice(const Matrix& g, const Matrix& p, double a, double b) {
    if (path.vertices.empty()) {
        path.vertices.push_back(current_location);
        visited_mask |= (1ULL << current_location);
    }

    if (current_location == end_location) { // Проверка на достижение конечной точки
        can_continue = false;
        return;
    }

    std::vector<std::size_t> neighbor_vertexes = getNeighborVertexes(g);

    // Если нет доступных соседей, завершаем путь
    if (neighbor_vertexes.empty()) {
        can_continue = false;
        return;
    }

    // Подсчет вероятности перехода муравья из одной вершины в другую
    std::vector<double> wish(neighbor_vertexes.size());
    double summary_wish = 0.0;

    for (std::size_t i = 0; i < neighbor_vertexes.size(); ++i) {
        std::size_t v = neighbor_vertexes[i];
        if (current_location >= g.GetRows() || v >= g.GetRows()) {
            throw std::out_of_range("Index is out of bounds in MakeChoice");
        }

        double t = p[current_location][v];
        double w = g[current_location][v];
        if (w == 0) {
            throw std::invalid_argument("Distance cannot be zero");
        }

        double n = 1 / w;
        wish[i] = std::pow(t, a) * std::pow(n, b);
        summary_wish += wish[i];
    }

    if (summary_wish <= 0.0) {
        can_continue = false;
        return;
    }

    // Подсчет вероятностей перехода
    std::vector<double> probability(neighbor_vertexes.size());
    for (std::size_t i = 0; i < neighbor_vertexes.size(); ++i) {
        probability[i] = wish[i] / summary_wish;
    }

    // Накопление вероятностей для бинарного поиска
    std::vector<double> choosing_probability(neighbor_vertexes.size());
    choosing_probability[0] = probability[0];
    for (std::size_t i = 1; i < neighbor_vertexes.size(); ++i) {
        choosing_probability[i] = choosing_probability[i - 1] + probability[i];
    }

    // Определение следующей вершины
    double choose = getRandomChoice();
    auto it = std::lower_bound(choosing_probability.begin(), choosing_probability.end(), choose);
    std::size_t next_vertex = neighbor_vertexes[it - choosing_probability.begin()];

    path.vertices.push_back(next_vertex);
    path.distance += g[current_location][next_vertex];
    visited_mask |= (1ULL << next_vertex);
    current_location = next_vertex;
}

aco_short_way::AntPath aco_short_way::AntColonyOptimization::SolveSalesmansProblem() {
    // Сколько итераций подряд можно без улучшений лучшего решения
    constexpr std::size_t kMaxIterationsWithoutImprovement = 200;

    const std::size_t n = graph_->GetRows();
    AntPath best_path;
    best_path.distance = std::numeric_limits<double>::infinity();

    std::size_t no_improve = 0;

    while (no_improve < kMaxIterationsWithoutImprovement) {
        // Локальная матрица прироста феромона за итерацию
        Matrix local_pheromone_update(n, n);  // по умолчанию заполнится нулями

        // Прогоняем всех муравьёв
        for (auto& ant : ants_) {
            // Сброс состояния муравья
            ant.visited_mask = 0;
            ant.can_continue = true;
            ant.current_location = start_location_;
            ant.path.vertices.clear();
            ant.path.distance = 0.0;

            // Двигаем муравья, пока он может идти
            while (ant.can_continue) {
                ant.MakeChoice(*graph_, *pheromone_, kAlpha_, kBeta_);
            }

            // Если муравей дошёл до цели — учитываем этот путь
            if (!ant.path.vertices.empty() &&
                ant.path.vertices.back() == end_location_ &&
                std::isfinite(ant.path.distance))
            {
                // Обновление лучшего решения
                if (ant.path.distance < best_path.distance) {
                    best_path = ant.path;
                    no_improve = 0;  // сброс стагнации
                }

                // Локальный депозит феромона по рёбрам маршрута
                // Чем короче путь — тем больше вклад
                const double deposit = kQ_ / std::max(1.0, ant.path.distance);
                for (std::size_t i = 0; i + 1 < ant.path.vertices.size(); ++i) {
                    const auto u = ant.path.vertices[i];
                    const auto v = ant.path.vertices[i + 1];
                    // Защита индексов (на всякий случай)
                    if (u < n && v < n) {
                        local_pheromone_update[u][v] += deposit;
                        local_pheromone_update[v][u] += deposit; // если граф неориентированный
                    }
                }
            }
        } // for ants

        // Глобальное обновление феромонов: испарение + депозит
        UpdateGlobalPheromone(local_pheromone_update);

        // Если в эту итерацию улучшения не было — считаем стагнацию
        ++no_improve;
    } // while

    return best_path;
}
