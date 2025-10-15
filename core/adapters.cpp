#include "registry.h"
#include <nlohmann/json.hpp>

// === Алгоритмы из твоего include/ ===
#include "../include/aco_tsp.h"
#include "../include/bbbc_tsp.h"
#include "../include/bbbc_hp.h"
#include "../include/aco_short.h"
#include "../include/Greedy_Adaptive_Method.h"
#include "../include/Dynamic_Grid_Method.h"
#include "../include/coa_tsp.h"        
#include "../src/algorithms/ss_hp.cpp"
#include "../src/algorithms/ss_tsp.cpp"

using namespace core;
using my_matrix::Matrix;

// ---------- ACO TSP ----------
static Result run_aco_tsp_impl(const Matrix& m, const Params& p) {
    std::size_t evaporation = p.value("evaporation", (std::size_t)1);
    std::size_t initial_pheromone = p.value("initial_pheromone", (std::size_t)1);
    std::size_t q_param = p.value("q", (std::size_t)5);
    std::size_t max_iterations = p.value("max_iter", (std::size_t)200);
    std::size_t max_no_improve = p.value("max_iter_no_improve", (std::size_t)48);
    std::size_t elite_ants = p.value("elite_ants", (std::size_t)5);

    aco_tsp::AntColonyOptimization aco(const_cast<Matrix&>(m),
        evaporation, initial_pheromone, q_param);
    auto best = aco.SolveSalesmansProblem(max_iterations, max_no_improve, elite_ants);

    Result r;
    r.algorithm = "aco_tsp"; r.problem = "tsp";
    r.cost = best.distance;
    for (auto v : best.vertices) r.path.push_back((int)v);
    r.meta = {
        {"evaporation", evaporation},
        {"initial_pheromone", initial_pheromone},
        {"q", q_param},
        {"max_iter", max_iterations},
        {"max_iter_no_improve", max_no_improve},
        {"elite_ants", elite_ants}
    };
    return r;
}

REGISTER_ALGO(aco_tsp, "tsp",
    nlohmann::json::parse(R"({
        "description": "Ant Colony Optimization for TSP",
        "defaults": {
            "evaporation": 1,
            "initial_pheromone": 1,
            "q": 5,
            "max_iter": 200,
            "max_iter_no_improve": 48,
            "elite_ants": 5
        }
    })"),
    run_aco_tsp_impl
);

// ---------- BBBC TSP ----------
static Result run_bbbc_tsp_impl(const Matrix& m, const Params& p) {
    int population = p.value("population", 200);
    int max_iters = p.value("max_iters", 2000);
    int elite = p.value("elite", 10);
    double mutation = p.value("mutation", 0.05);
    unsigned seed = p.value("seed", 42u);
    size_t max_neighbors = p.value("max_neighbors", (size_t)20);

    bbbc::BigBangBigCrunch bbbc(const_cast<Matrix&>(m),
        population, max_iters, elite, mutation, seed, max_neighbors);
    bbbc.solve();

    Result r;
    r.algorithm = "bbbc_tsp"; r.problem = "tsp";
    r.cost = bbbc.getBestDistance();
    auto& route = bbbc.getBestRoute();
    r.path.assign(route.begin(), route.end());
    r.meta = {
        {"population", population},
        {"max_iters", max_iters},
        {"elite", elite},
        {"mutation", mutation},
        {"seed", seed},
        {"max_neighbors", max_neighbors}
    };
    return r;
}

REGISTER_ALGO(bbbc_tsp, "tsp",
    nlohmann::json::parse(R"({
        "description": "Big Bang Big Crunch for TSP",
        "defaults": {
            "population": 200, "max_iters": 2000, "elite": 10,
            "mutation": 0.05, "seed": 42, "max_neighbors": 20
        }
    })"),
    run_bbbc_tsp_impl
);

// ---------- BBBC HP (Hamiltonian Path со старт/финиш) ----------
static Result run_bbbc_hp_impl(const Matrix& m, const Params& p) {
    int population = p.value("population", 200);
    int max_iters = p.value("max_iters", 2000);
    int elite = p.value("elite", 10);
    double mutation = p.value("mutation", 0.05);
    int start = p.value("start", 0);     // 0-based
    int end = p.value("end", (int)m.GetRows() - 1);
    unsigned seed = p.value("seed", 42u);
    size_t max_neighbors = p.value("max_neighbors", (size_t)20);

    bbbc_hp::BigBangBigCrunch_HP hp(const_cast<Matrix&>(m),
        population, max_iters, elite, mutation, start, end, seed, max_neighbors);
    hp.solve();

    Result r;
    r.algorithm = "bbbc_hp"; r.problem = "hp";
    r.cost = hp.getBestDistance();
    auto& route = hp.getBestRoute();
    r.path.assign(route.begin(), route.end());
    r.meta = {
        {"population", population},
        {"max_iters", max_iters},
        {"elite", elite},
        {"mutation", mutation},
        {"start", start},
        {"end", end},
        {"seed", seed},
        {"max_neighbors", max_neighbors}
    };
    return r;
}

REGISTER_ALGO(bbbc_hp, "hp",
    nlohmann::json::parse(R"({
        "description": "Big Bang Big Crunch for Hamiltonian Path [start/end]",
        "defaults": {
    "population":200, "max_iters" : 2000, "elite" : 10,
        "mutation" : 0.05, "start" : 0, "end" : -1, "seed" : 42,
        "max_neighbors" : 20
}
    })"),
    run_bbbc_hp_impl
        );

// ---------- ACO Shortest Path (между start и end) ----------
static Result run_aco_short_impl(const Matrix& m, const Params& p) {
    std::size_t start = p.value("start", (std::size_t)0);
    std::size_t end = p.value("end", (std::size_t)std::max(1u, m.GetRows()) - 1);

    aco_short_way::AntColonyOptimization aco(const_cast<Matrix&>(m), start, end);
    auto best = aco.SolveShortestPath();

    Result r;
    r.algorithm = "aco_short"; r.problem = "shortest_path";
    r.cost = best.distance;
    for (auto v : best.vertices) r.path.push_back((int)v);
    r.meta = { {"start", start}, {"end", end} };
    return r;
}

REGISTER_ALGO(aco_short, "shortest_path",
    nlohmann::json::parse(R"({
        "description": "Ant Colony Optimization for shortest path",
        "defaults": { "start": 0, "end" : -1 }
    })"),
    run_aco_short_impl
);

// ---------- Greedy Adaptive (TSP) ----------
static Result run_greedy_adaptive_tsp_impl(const Matrix& m, const Params& p) {
    size_t iters = p.value("iterations", (size_t)500);
    greedy_adaptive::GreedyAdaptiveSolver solver(m, iters);
    auto sol = solver.solve();

    Result r;
    r.algorithm = "greedy_tsp"; r.problem = "tsp";
    r.cost = sol.distance;
    for (auto v : sol.path) r.path.push_back((int)v);
    r.meta = { {"iterations", iters} };
    return r;
}

REGISTER_ALGO(greedy_tsp, "tsp",
    nlohmann::json::parse(R"({
        "description": "Greedy Adaptive with 2-opt for TSP",
        "defaults": { "iterations": 500 }
    })"),
    run_greedy_adaptive_tsp_impl
);

// ---------- Dynamic Grid (TSP) ----------
static Result run_grid_tsp_impl(const Matrix& m, const Params& p) {
    dynamic_grid::GridParams gp;
    gp.initial_population = p.value("initial_population", (size_t)200);
    gp.max_iterations = p.value("max_iterations", (size_t)1000);
    gp.max_local_iter = p.value("max_local_iter", (size_t)200);
    gp.local_expansion_factor = p.value("local_expansion_factor", 0.5);  // доля от текущего размера популяции
    gp.global_expansion_factor = p.value("global_expansion_factor", 0.2);
    gp.use_3opt = p.value("use_3opt", false);
    gp.parallel = p.value("parallel", false);
    gp.verbose = p.value("verbose", false);
    gp.print_interval = p.value("print_interval", (size_t)100);

    dynamic_grid::DynamicGridSolver solver(m, gp);
    auto sol = solver.solve();

    Result r;
    r.algorithm = "grid_tsp"; r.problem = "tsp";
    r.cost = sol.distance;
    for (auto v : sol.path) r.path.push_back((int)v);
    r.meta = {
        {"initial_population", gp.initial_population},
        {"max_iterations", gp.max_iterations},
        {"max_local_iter", gp.max_local_iter},
        {"local_expansion_factor", gp.local_expansion_factor},
        {"global_expansion_factor", gp.global_expansion_factor},
        {"use_3opt", gp.use_3opt},
        {"parallel", gp.parallel},
        {"verbose", gp.verbose},
        {"print_interval", gp.print_interval}
    };
    return r;
}

REGISTER_ALGO(grid_tsp, "tsp",
    nlohmann::json::parse(R"({
        "description": "Dynamic Grid metaheuristic for TSP",
        "defaults": {
            "initial_population":200, "max_iterations":1000, "max_local_iter":200,
            "local_expansion_factor":0.5, "global_expansion_factor":0.2,
            "use_3opt": false, "parallel": false, "verbose": false, "print_interval":100
        }
    })"),
    run_grid_tsp_impl
);

// ---------- Cuckoo / COA (TSP) ----------
static Result run_coa_tsp_impl(const Matrix& m, const Params& p) {
    size_t population = p.value("population", (size_t)100);
    double pa = p.value("pa", 0.25);          // discovery probability
    double step = p.value("step_size", 1.0);    // шаг Леви
    int max_iters = p.value("max_iters", 1000);

    coa_tsp::CuckooSearch cs(const_cast<Matrix&>(m), population, pa, step, max_iters);
    auto path = cs.solveTSP();

    // считаем длину по матрице (на всякий)
    double cost = 0.0;
    if (!path.empty()) {
        for (size_t i = 0; i + 1 < path.size(); ++i) cost += m[path[i]][path[i + 1]];
        cost += m[path.back()][path.front()];
    }

    Result r;
    r.algorithm = "coa_tsp"; r.problem = "tsp";
    r.cost = cost;
    for (auto v : path) r.path.push_back((int)v);
    r.meta = {
        {"population", population}, {"pa", pa},
        {"step_size", step}, {"max_iters", max_iters}
    };
    return r;
}

REGISTER_ALGO(coa_tsp, "tsp",
    nlohmann::json::parse(R"({
        "description": "Cuckoo Search for TSP",
        "defaults": { "population":100, "pa":0.25, "step_size":1.0, "max_iters":1000 }
    })"),
    run_coa_tsp_impl
);

// ---------- SS TSP (reference set; 1-based внутри, аккуратно) ----------
static Result run_ss_tsp_impl(const Matrix& m, const Params& p) {
    bool debug = p.value("debug", false);
    unsigned iterations = p.value("iterations", 150u);
    unsigned reference_size = p.value("reference_size", 10u);
    unsigned elite_size = p.value("elite_size", 2u);
    double reverse_prob = p.value("reverse_prob", 0.5);
    double scramble_prob = p.value("scramble_prob", 0.3);

    // ВНИМАНИЕ: реализации ss_* оперируют вершинами 1..N и кладут последний=первый.
    // Они принимают const Matrix*.
    auto* best = StartHamiltonCycle(
        debug, &m, iterations, reference_size, elite_size, reverse_prob, scramble_prob);

    Result r;
    r.algorithm = "ss_tsp"; r.problem = "tsp";
    r.cost = best->distance;

    // Перегоняем CustomVector<double> (1..N) -> 0..N-1, цикл оставляем
    // (ожидаем, что размер = N+1, где последний = первый)
    if (best->tour) {
        for (int i = 0; i < best->tour->GetSize(); ++i) {
            int v1 = (int)(*(best->tour))[i] - 1;
            r.path.push_back(v1);
        }
    }
    r.meta = {
        {"iterations", iterations}, {"reference_size", reference_size},
        {"elite_size", elite_size}, {"reverse_prob", reverse_prob},
        {"scramble_prob", scramble_prob}
    };

    // освобождение best/его внутренних — в их коде это делаетсянаружи;
    // если нужно, добавь безопасное удаление, когда будешь рефакторить ss_*.

    return r;
}

REGISTER_ALGO(ss_tsp, "tsp",
    nlohmann::json::parse(R"({
        "description": "SS (reference set) for TSP [1-based internally]",
        "defaults": {
            "iterations":150, "reference_size" : 10, "elite_size" : 2,
            "reverse_prob" : 0.5, "scramble_prob" : 0.3, "debug" : false
        }
    })"),
    run_ss_tsp_impl
);

// ---------- SS HP (1-based внутри; start/end внешне в 0-based) ----------
static Result run_ss_hp_impl(const Matrix& m, const Params& p) {
    bool debug = p.value("debug", false);
    unsigned iterations = p.value("iterations", 150u);
    unsigned reference_size = p.value("reference_size", 10u);
    double reverse_prob = p.value("reverse_prob", 0.5);
    double scramble_prob = p.value("scramble_prob", 0.3);
    int start0 = p.value("start", 0);
    int end0 = p.value("end", (int)m.GetRows() - 1);

    // SS HP ждёт 1..N, поэтому переводим:
    int start1 = start0 + 1;
    int end1 = end0 + 1;

    auto* best = StartHamiltonPath(
        debug, &m, iterations, reference_size, reverse_prob, scramble_prob, start1, end1);

    Result r;
    r.algorithm = "ss_hp"; r.problem = "hp";
    r.cost = best->distance;

    if (best->tour) {
        for (int i = 0; i < best->tour->GetSize(); ++i) {
            r.path.push_back((int)(*(best->tour))[i] - 1); // назад в 0-based
        }
    }
    r.meta = {
        {"iterations", iterations}, {"reference_size", reference_size},
        {"reverse_prob", reverse_prob}, {"scramble_prob", scramble_prob},
        {"start", start0}, {"end", end0}, {"debug", debug}
    };

    return r;
}

REGISTER_ALGO(ss_hp, "hp",
    nlohmann::json::parse(R"({
        "description": "SS for Hamiltonian Path [1-based internally]",
        "defaults": {
            "iterations":150, "reference_size" : 10,
            "reverse_prob" : 0.5, "scramble_prob" : 0.3,
            "start" : 0, "end" : -1, "debug" : false
        }
    })"),
    run_ss_hp_impl
);

extern "C" void force_link_adapters() {}
