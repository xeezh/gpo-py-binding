#include <iostream>
#include <list>
#include <algorithm>
#include "random"

#include "../classes/CustomVector.h"
#include "../../include/Solution.h"

using namespace structures;
using namespace std;

CustomVector<double> *CreateTourHamiltonCycle(const Matrix *distance_matrix) {
    unsigned int city_num = distance_matrix->GetRows();
    CustomVector<double> *out_tour = new CustomVector<double>( city_num + 1);

    for (int i = 0; i < city_num; ++i) {
        (*out_tour)[i] = i + 1;
    }
    (*out_tour)[city_num] = (*out_tour)[0];
    return out_tour;
}

Solution *RandomizeTourHamiltonCycle(const Matrix *distance_matrix, CustomVector<double> *input_tour) {
    CustomVector<double> *out_tour = new CustomVector<double>(*input_tour);
    out_tour->Shuffle(0, out_tour->GetSize() - 1);
    (*out_tour)[out_tour->GetSize() - 1] = (*out_tour)[0];
    Solution *res = new Solution();
    res->tour = out_tour;
    res->distance = TourDistanceCalc(distance_matrix, out_tour);
    return res;
}

void LocalSearchHamiltonCycle(const Matrix *distance_matrix, Solution *sol) {
    unsigned int sol_len = sol->tour->GetSize();

    Solution best_route = Solution(*sol);
    Solution seed = Solution(*sol);

    for (unsigned int i = 0; i < sol_len - 2; ++i) {
        for (unsigned int j = i + 1; j < sol_len - 1; ++j) {
            best_route.tour->Reverse(i, j);
            (*(best_route.tour))[sol_len - 1] = (*(best_route.tour))[0];
            best_route.distance = TourDistanceCalc(distance_matrix, best_route.tour);
            if (best_route.distance < sol->distance) {
                *sol = best_route;
            }
            best_route = seed;
        }
    }
}

Solution *CrossoverHamiltonCycle(const Matrix *distance_matrix,
                                 const std::pair<Solution *, Solution *> &parents,
                                 double reverse_prob = 0.5,
                                 double scramble_prob = 0.3) {
    std::random_device generator;
    std::uniform_real_distribution<double> double_distribution(0, 1);

    CustomVector<double> *parent1 = new CustomVector<double>(*parents.first->tour);
    CustomVector<double> *parent2 = new CustomVector<double>(*parents.second->tour);

    unsigned int i = std::rand() % (parent1->GetSize() - 2);
    unsigned int j = std::rand() % (parent1->GetSize() - 2);
    if (i > j) std::swap(i, j);

    double chance = double_distribution(generator);
    if (chance < reverse_prob) { // Прокнул реверс
        parent1->Reverse(i, j);
    }

    unsigned int insert_index = 0;
    bool can_insert;
    for (unsigned int k = 0; k < parent2->GetSize() - 1; ++k) {
        can_insert = true;
        for (unsigned int l = i; l <= j; ++l) {
            if ((*parent2)[k] == (*parent1)[l]) {
                can_insert = false;
                break;
            }
        }
        if (can_insert) {
            (*parent2)[insert_index] = (*parent2)[k];
            insert_index++;
        }
    }

    chance = double_distribution(generator);
    if (chance < scramble_prob && insert_index != 0) {
        parent2->Shuffle(0, insert_index);
    }

    insert_index = 0;
    for (unsigned int k = 0; k < i; ++k) {
        (*parent1)[k] = (*parent2)[insert_index];
        insert_index++;
    }
    for (unsigned int k = j + 1; k < parent1->GetSize() - 1; ++k) {
        (*parent1)[k] = (*parent2)[insert_index];
        insert_index++;
    }

    (*parent1)[parent1->GetSize() - 1] = (*parent1)[0];
    delete parent2;
    Solution *sol = new Solution();
    sol->tour = parent1;
    sol->distance = TourDistanceCalc(distance_matrix, parent1);
    return sol;
}

Solution *StartHamiltonCycle(
        bool input_debug,
        const Matrix *dist_matrix,
        unsigned int _iterations = 150,
        unsigned int _reference_size = 10,
        unsigned int _elite_size = 2,
        double _reverse_prob = 0.5,
        double _scramble_prob = 0.3) {
    unsigned int iterations = _iterations;
    unsigned int reference_size = _reference_size;
    unsigned int elite_size = _elite_size;
    double reverse_prob = _reverse_prob;
    double scramble_prob = _scramble_prob;
    bool debug = input_debug;

    unsigned int city_num = dist_matrix->GetRows();
    CustomVector<double> *city_tour = CreateTourHamiltonCycle(dist_matrix);
    double dist = TourDistanceCalc(dist_matrix, city_tour);

    unsigned int iter = 0;
    Solution *best = new Solution(city_tour, dist);
    list<Solution *> reference_list;
    list<Solution *> elite_list;

    for (int i = 0; i < reference_size; ++i) {
        auto tmp = RandomizeTourHamiltonCycle(dist_matrix, city_tour);
        reference_list.push_back(tmp);
    }

    while (iter < iterations) {
        if (debug && iter % 15 == 0) {
            cout << "Iteration: " << iter + 1 << " Best sol: " << best->distance << endl;
        }

        iter++;
        list<Solution *> candidate_list;

        reference_list.sort(CompareSolutions);
        elite_list.clear();
        auto it = reference_list.begin();
        for (unsigned int i = 0; i < elite_size && it != reference_list.end(); ++i, ++it) {
            elite_list.push_back(new Solution(**it));
        }

        auto elite_it = elite_list.begin();
        auto ref_it = reference_list.begin();

        while (elite_it != elite_list.end() && ref_it != reference_list.end()) {
            while (std::find(elite_list.begin(), elite_list.end(), *ref_it) != elite_list.end()) {
                ++ref_it;
                if (ref_it == reference_list.end()) break;
            }
            if (ref_it == reference_list.end()) break;

            auto child = CrossoverHamiltonCycle(dist_matrix, { *elite_it, *ref_it }, reverse_prob, scramble_prob);
            candidate_list.push_back(child);
            ++elite_it;
            ++ref_it;
        }

        for (auto &candidate : candidate_list) {
            LocalSearchHamiltonCycle(dist_matrix, candidate);
        }

        for (auto candidate : candidate_list) {
            reference_list.push_back(candidate);
        }
        reference_list.sort(CompareSolutions);
        while (reference_list.size() > reference_size) {
            delete reference_list.back();
            reference_list.pop_back();
        }

        for (auto sol : reference_list) {
            if (sol->distance < best->distance) {
                *best = *sol;
            }
        }

        candidate_list.clear();
    }

    while (!reference_list.empty()) {
        delete reference_list.back();
        reference_list.pop_back();
    }
    while (!elite_list.empty()) {
        delete elite_list.back();
        elite_list.pop_back();
    }

    delete city_tour;
    return best;
}
