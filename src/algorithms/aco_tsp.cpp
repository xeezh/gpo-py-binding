#include <iostream>
#include <atomic>
#include <random>
#include <thread>
#include "../../include/aco_tsp.h"
#include <chrono>

using namespace aco_tsp;

/* Конструктор колонии */
aco_tsp::AntColonyOptimization::AntColonyOptimization(Matrix& graph, std::size_t evaporation,
                                                      std::size_t initial_pheromone,
                                                      std::size_t q) : kQ_(q * graph.GetRows()), graph_(&graph),
                                                      pheromone_(nullptr), kEvaporation_(evaporation), kPheromone0_(initial_pheromone) {
    const std::size_t kVertexesCount = graph_->GetRows();

    std::vector<double> array(kVertexesCount * kVertexesCount, kPheromone0_);
    pheromone_ = std::make_unique<Matrix>(kVertexesCount, kVertexesCount, array.data());
    std::cout << "Pheromone Matrix initialized with size: " << kVertexesCount << "x" << kVertexesCount << "\n";

    CreateAnts(); // создание муравьев
    std::cout << "Created " << kVertexesCount << " ants.\n";
}

/* Инициализация муравьев */
void aco_tsp::AntColonyOptimization::CreateAnts() {
    const auto kAntsCount = graph_->GetRows();
    ants_.resize(kAntsCount);

    for (std::size_t i = 0; i != kAntsCount; ++i)
        ants_[i] = Ant(i);
}

/* Локальное обновление феромонов */
void aco_tsp::AntColonyOptimization::UpdateLocalPheromone(std::size_t from, std::size_t to) {
    if (!pheromone_ || pheromone_->GetRows() == 0) {
        throw std::runtime_error("Pheromone Matrix is not initialized");
    }
    if (from >= pheromone_->GetRows() || to >= pheromone_->GetRows()) {
        std::cerr << "Error: Index out of bounds in UpdateLocalPheromone: from=" << from << ", to=" << to << "\n";
        throw std::out_of_range("Index out of bounds in UpdateLocalPheromone");
    }

    (*pheromone_)[from][to] = (1 - kEvaporation_) * (*pheromone_)[from][to] + 0.01;
    //std::cout << "kEvaporation_: " << kEvaporation_ << "\n";
    (*pheromone_)[to][from] = (*pheromone_)[from][to]; // Симметрия
}

/* Глобальное обновление феромонов */
void aco_tsp::AntColonyOptimization::UpdateGlobalPheromone(const Matrix& lpu, const std::vector<AntPath>& elite_paths) const {
    for (std::size_t from = 0; from < lpu.GetRows(); ++from) {
        for (std::size_t to = 0; to < lpu.GetRows(); ++to) {
            (*pheromone_)[from][to] = (1 - kEvaporation_) * (*pheromone_)[from][to] + lpu[from][to];
        }
    }
    // Усиление феромонов для элитных путей
    for (const auto& path : elite_paths) {
        if (path.vertices.size() < 2) continue; // Пропускаем некорректные пути
        double elite_pheromone = kQ_ / path.distance;
        for (std::size_t i = 0; i < path.vertices.size() - 1; ++i) {
            std::size_t from = path.vertices[i];
            std::size_t to = path.vertices[i + 1];
            if (from < pheromone_->GetRows() && to < pheromone_->GetRows()) {
                (*pheromone_)[from][to] += elite_pheromone;
                (*pheromone_)[to][from] = (*pheromone_)[from][to]; // Симметрия
            }
        }
    }
}

/* Динамическое регулирование параметров */
void aco_tsp::AntColonyOptimization::AdjustParameters(std::size_t iteration, std::size_t max_iterations) {
    double progress = static_cast<double>(iteration) / max_iterations;
    kAlpha_ = 1.5 + progress * 1.0; // Увеличиваем от 1.5 до 2.5
    kBeta_ = 2.5 - progress * 1.0; // Уменьшаем от 2.5 до 1.5
}

/* Получение случайного числа в пределах (0:1) */
double aco_tsp::Ant::getRandomChoice() {
    std::uniform_real_distribution<> dist(0.0, 1.0);
    std::default_random_engine re(std::chrono::system_clock::now().time_since_epoch().count());
    return dist(re);
}

/* Поиск соседних вершин */
std::vector<std::size_t> aco_tsp::Ant::getNeighborVertexes(const Matrix& g) const {
    std::vector<std::size_t> neighbors;
    for (std::size_t i = 0; i < g.GetRows(); ++i) {
        if (g[current_location][i] > 0 && !(visited_mask & (1ULL << i))) {
            neighbors.push_back(i);
        }
    }
    return neighbors;
}

/* Функция выбора вершины муравьем */
void aco_tsp::Ant::MakeChoice(const Matrix& g, const Matrix& p, double a, double b) {
    if (current_location >= g.GetRows()) {
        std::cerr << "Error: current_location out of bounds: " << current_location << "\n";
        throw std::out_of_range("Invalid current_location in MakeChoice");
    }

    // Получаем соседей
    std::vector<std::size_t> neighbor_vertexes = getNeighborVertexes(g);

    if (neighbor_vertexes.empty()) {
        can_continue = false;
        return;
    }

    // Если все вершины посещены, вернуться на начальную точку
    if (visited_mask == (1ULL << g.GetRows()) - 1) {
        // Находим начальную вершину
        for (std::size_t i = 0; i < g.GetRows(); ++i) {
            if (!(visited_mask & (1ULL << i))) {
                current_location = i;
                visited_mask |= (1ULL << i);
                break;
            }
        }
        return;
    }

    // Рассчитываем вероятности перехода
    std::vector<double> wish(neighbor_vertexes.size());
    double summary_wish = 0.0;

    for (std::size_t i = 0; i < neighbor_vertexes.size(); ++i) {
        std::size_t v = neighbor_vertexes[i];
        if (current_location >= g.GetRows() || v >= g.GetRows()) {
            std::cerr << "Error: Neighbor vertex index out of bounds: v=" << v << "\n";
            throw std::out_of_range("Invalid neighbor vertex index in MakeChoice");
        }

        double t = p[current_location][v];
        double w = g[current_location][v];
        if (w == 0) {
            std::cerr << "Error: Zero weight for edge (" << current_location << ", " << v << ")\n";
            throw std::invalid_argument("Zero weight edge in MakeChoice");
        }

        double n = 1.0 / w;
        wish[i] = std::pow(t, a) * std::pow(n, b);
        summary_wish += wish[i];
    }

    if (summary_wish == 0) {
        std::cerr << "Error: Summary wish is zero. No valid moves possible.\n";
        throw std::runtime_error("No valid moves in MakeChoice");
    }

    // Рассчитываем вероятности для выбора вершины
    std::vector<double> probability(neighbor_vertexes.size());
    for (std::size_t i = 0; i < neighbor_vertexes.size(); ++i) {
        probability[i] = wish[i] / summary_wish;
    }

    std::vector<double> choosing_probability(neighbor_vertexes.size());
    choosing_probability[0] = probability[0];
    for (std::size_t i = 1; i < neighbor_vertexes.size(); ++i) {
        choosing_probability[i] = choosing_probability[i - 1] + probability[i];
    }

    double choose = getRandomChoice();
    auto it = std::lower_bound(choosing_probability.begin(), choosing_probability.end(), choose);
    if (it == choosing_probability.end()) {
        std::cerr << "Error: Random choice out of bounds: " << choose << "\n";
        throw std::out_of_range("Random choice index out of bounds in MakeChoice");
    }

    // Выбираем следующую вершину
    current_location = neighbor_vertexes[it - choosing_probability.begin()];
    visited_mask |= (1ULL << current_location); // Отмечаем посещение вершины
}

/* Сброс феромонов */
void aco_tsp::AntColonyOptimization::ResetPheromones() {
    const std::size_t kVertexesCount = pheromone_->GetRows();
    for (std::size_t i = 0; i < kVertexesCount; ++i) {
        for (std::size_t j = 0; j < kVertexesCount; ++j) {
            (*pheromone_)[i][j] = kPheromone0_;
        }
    }
}

/* Задача коммивояжера */
aco_tsp::AntPath aco_tsp::AntColonyOptimization::SolveSalesmansProblem(
    std::size_t max_iterations,
    std::size_t max_iterations_without_improvement,
    std::size_t elite_ants_count) {
    const std::size_t kVertexesCount = graph_->GetRows();

    AntPath best_path;
    best_path.distance = 1e9;

    std::vector<AntPath> elite_paths(elite_ants_count, AntPath{std::vector<std::size_t>(), 1e9});

    if (!graph_ || graph_->GetRows() == 0 || graph_->GetCols() != graph_->GetRows()) {
        throw std::runtime_error("Invalid input graph: Matrix is null or not square");
    }

    if (!pheromone_ || pheromone_->GetRows() != graph_->GetRows()) {
        throw std::runtime_error("Invalid pheromone Matrix: not initialized or size mismatch");
    }

    std::atomic<std::size_t> counter(0);
    for (std::size_t iteration = 0; iteration < max_iterations; ++iteration) {
        Matrix local_pheromone_update(kVertexesCount, kVertexesCount);

        // Обход всех муравьев
        for (auto& ant : ants_) {
            std::vector<std::size_t> path;
            double path_distance = 0.0;

            path.push_back(ant.start_location);
            ant.visited_mask = (1ULL << ant.start_location); // Начинаем с начальной вершины

            ant.current_location = ant.start_location;
            ant.can_continue = true;

            // Строим путь муравья
            while (ant.can_continue) {
                std::size_t previous_location = ant.current_location;
                ant.MakeChoice(*graph_, *pheromone_, kAlpha_, kBeta_);

                if (ant.can_continue) {
                    path.push_back(ant.current_location);
                    path_distance += (*graph_)[previous_location][ant.current_location];
                    UpdateLocalPheromone(previous_location, ant.current_location);
                }
            }

            // Возвращаемся на начальную вершину
            path.push_back(ant.start_location);
            path_distance += (*graph_)[ant.current_location][ant.start_location];

            // Проверка пути муравья
            if (path.size() == kVertexesCount + 1) {
                if (best_path.distance > path_distance) {
                    best_path = AntPath(path, path_distance);
                    counter.store(0); // Сброс счетчика
                    //std::cout << "Current best distance: " << best_path.distance << "\n";
                }

                // Обновление элитных путей
                for (auto& elite_path : elite_paths) {
                    if (elite_path.distance > path_distance) {
                        elite_path = AntPath(path, path_distance);
                        break;
                    }
                }
            }
        }

        // Глобальное обновление феромонов с использованием элитных путей
        UpdateGlobalPheromone(local_pheromone_update, elite_paths);

        // Динамическое регулирование параметров
        AdjustParameters(iteration, max_iterations);

        // Проверка на остановку или сброс
        if (++counter >= max_iterations_without_improvement) {
            std::cout << "Resetting pheromones due to lack of improvement\n";
            ResetPheromones(); // Сброс феромонов
            counter.store(0);  // Сброс счетчика
        }
    }
    return best_path;
}

