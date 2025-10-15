#include <iostream>
#include <unordered_set>
#include "../../include/Solution.h"
#include "../../include/bbbc_tsp.h"

using std::cout, std::endl;
using namespace bbbc;
using namespace structures;

// Конструктор
BigBangBigCrunch::BigBangBigCrunch(Matrix& graph, int populationSize, int maxIterations, int eliteSize, double mutationRate, unsigned int seed, size_t maxNeighbors)
        : graph_(graph),
          Population_Size(populationSize),
          Max_Iterations(maxIterations),
          Elite_Size(eliteSize),
          Mutation_Rate(mutationRate),
          bestDistance(std::numeric_limits<double>::max()),
          rng(seed),
          Max_Neighbors(maxNeighbors)
{
    initializePopulation();

    elitePool.resize(Elite_Size);
    eliteDistances.resize(Elite_Size, std::numeric_limits<double>::max());
}

double BigBangBigCrunch::getBestDistance() const {
    return bestDistance;
}

const std::vector<int> &BigBangBigCrunch::getBestRoute() const {
    return bestRoute;
}

double BigBangBigCrunch::Tour_Distance_Calc(const std::vector<int>& city_tour) {
    double totalDistance = 0.0;

    for (size_t i = 1; i < city_tour.size(); ++i) {
        totalDistance += graph_[city_tour[i - 1]][city_tour[i]];
    }
    totalDistance += graph_[city_tour.back()][city_tour.front()];


    return totalDistance;
}


void BigBangBigCrunch::initializePopulation() {
    population.clear();

    // Создаём базовый маршрут 0,1,2,...,N-1
    std::vector<int> baseTour(graph_.GetRows());
    for (size_t i = 0; i < baseTour.size(); ++i) {
        baseTour[i] = static_cast<int>(i);
    }

    for (int i = 0; i < Population_Size; ++i) {
        // Копируем базовый тур и перемешиваем случайно
        std::vector<int> tour = baseTour;
        std::shuffle(tour.begin(), tour.end(), rng);

        population.push_back(tour);

        double dist = Tour_Distance_Calc(tour);
        if (dist < bestDistance) {
            bestDistance = dist;
            bestRoute = tour;
        }
    }
}


// Создание соседей для текущего маршрута
std::vector<std::vector<int>> BigBangBigCrunch::createNeighbors(const std::vector<int>& route ) const {
    std::vector<std::vector<int>> neighbors;

    for (size_t i = 1; i < route.size() - 2 && neighbors.size() < Max_Neighbors; ++i) {
        for (size_t j = i + 1; j < route.size() - 1 && neighbors.size() < Max_Neighbors; ++j) {
            std::vector<int> newRoute(route);
            std::reverse(newRoute.begin() + i, newRoute.begin() + j + 1);
            neighbors.push_back(newRoute);
        }
    }
    return neighbors;
}


double BigBangBigCrunch::localSearch(std::vector<int>& route, size_t max_neighbors) {
    bool improved = true;
    double currentDistance = Tour_Distance_Calc(route);

    while (improved) {
        improved = false;
        auto neighbors = createNeighbors(route);

        double bestNeighborDist = currentDistance;
        std::vector<int> bestNeighbor = route;

        for (const auto& neighbor : neighbors) {
            double dist = Tour_Distance_Calc(neighbor);

            if (dist < bestNeighborDist) {
                bestNeighborDist = dist;
                bestNeighbor = neighbor;
            }
        }

        if (bestNeighborDist < currentDistance) {
            route = std::move(bestNeighbor);
            currentDistance = bestNeighborDist;
            improved = true;
        }
    }

    return currentDistance; // ⬅ теперь возвращает расстояние
}



// Обновление элитного пула
void BigBangBigCrunch::updateElitePool(const std::vector<int>& newRoute, double newDistance) {
    // Если лучше текущего лучшего — обновляем лучший результат
    if (newDistance < bestDistance) {
        bestDistance = newDistance;
        bestRoute = newRoute;
    }

    // Добавляем в элиту если лучше худшего в пуле
    auto maxIt = std::max_element(eliteDistances.begin(), eliteDistances.end());
    if (*maxIt > newDistance) {
        int idx = std::distance(eliteDistances.begin(), maxIt);
        eliteDistances[idx] = newDistance;
        elitePool[idx] = newRoute;
    }
}

void BigBangBigCrunch::removeWorstSolutions() {
    std::sort(population.begin(), population.end(), [&](const std::vector<int>& a, const std::vector<int>& b){
        return Tour_Distance_Calc(a) < Tour_Distance_Calc(b);
    });

    if (population.size() > Elite_Size) {
        population.resize(Elite_Size);
    }
}

bool BigBangBigCrunch::isValidTour(const std::vector<int>& tour) {
    if (tour.size() != graph_.GetRows()) return false;

    std::vector<bool> visited(tour.size(), false);
    for (int city : tour) {
        if (city < 0 || city >= static_cast<int>(tour.size()) || visited[city]) {
            return false; // город вне диапазона или уже посещён
        }
        visited[city] = true;
    }

    return true; // все города встречаются ровно один раз

}

// Метод для запуска алгоритма
void BigBangBigCrunch::solve() {
    for (int iteration = 0; iteration < Max_Iterations; ++iteration) {
        double totalDistance = 0.0;
        double worstDistance = std::numeric_limits<double>::lowest();

        // Фаза «сжатия»: улучшаем каждую особь локальным поиском
        for (auto& route : population) {
            double dist = localSearch(route, Max_Neighbors); // ← возвращает уже посчитанную дистанцию
            updateElitePool(route, dist);

            if (dist < bestDistance) {
                bestDistance = dist;
                bestRoute = route;
            }

            totalDistance += dist;
            if (dist > worstDistance) worstDistance = dist;
        }

        // Фаза «взрыва»: генерируем новых соседей из элитных решений
        std::vector<std::vector<int>> newPopulation;
        for (const auto& eliteRoute : elitePool) {
            auto neighbors = createNeighbors(eliteRoute);
            for (auto& neighbor : neighbors) {
                // Возможна мутация с вероятностью Mutation_Rate
                if (std::uniform_real_distribution<double>(0.0, 1.0)(rng) < Mutation_Rate) {
                    // Пример простой мутации — перестановка двух городов
                    size_t idx1 = std::uniform_int_distribution<size_t>(0, neighbor.size()-1)(rng);
                    size_t idx2 = std::uniform_int_distribution<size_t>(0, neighbor.size()-1)(rng);
                    std::swap(neighbor[idx1], neighbor[idx2]);
                }


                double dist =  localSearch(neighbor, Max_Neighbors);


                if (dist < bestDistance) {
                    bestDistance = dist;
                    bestRoute = neighbor;
                }

                newPopulation.push_back(neighbor);
            }
        }
        if (!isValidTour(bestRoute)) {
            std::cerr << "Invalid tour detected after mutation!" << std::endl;
        }

        // Если новая популяция пуста (например, на первых итерациях) — сохраняем старую
        if (!newPopulation.empty()) {
            population = std::move(newPopulation);
        }

        removeWorstSolutions();

        double avgDistance = totalDistance / population.size();

        if (iteration % 100 == 0 || iteration == Max_Iterations - 1) {
            cout << "Iteration " << iteration
                 << ": Best distance = " << bestDistance
                 << ", Avg distance = " << avgDistance
                 << ", Worst distance = " << worstDistance
                 << endl;
        }
    }
}

// Опциональный метод для вывода лучшего маршрута
void BigBangBigCrunch::printBestRoute() const {
    cout << "Best route: ";
    for (int city : bestRoute) {
        cout << city << " ";
    }
    cout << "\nDistance: " << bestDistance << endl;
}

// Деструктор — теперь ничего не удаляем, graph_ - ссылка
BigBangBigCrunch::~BigBangBigCrunch() {
    population.clear();
    elitePool.clear();
    eliteDistances.clear();
}