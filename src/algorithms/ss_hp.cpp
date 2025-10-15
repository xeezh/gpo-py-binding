#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <list>
#include <cmath>
#include <algorithm>
#include "random"
#include <ctime>
#include <limits>
#include "../../include/Matrix.h"
#include "../../include/problem_parser.h"

#include "../classes/CustomVector.h"
#include "../../include/Solution.h"

using namespace structures;
using namespace std;

CustomVector<double> *CreateTourHamiltonPath(const Matrix *distance_matrix, int start, int end) {
    /* start и end не должны отсчитываться от нуля */
    unsigned int city_num = distance_matrix->GetRows();
    CustomVector<double> *out_tour = new CustomVector<double>(city_num);
    (*out_tour)[0] = start;
    (*out_tour)[city_num-1] = end;
    int insert_index = 1;

    for (int i = 1; i < start; ++i) {
        (*out_tour)[insert_index++] = i;
    }
    for (int i = start+1; i < end; ++i) {
        (*out_tour)[insert_index++] = i;
    }
    for (int i = end+1; i <= city_num; ++i) {
        (*out_tour)[insert_index++] = i;
    }
    return out_tour;
}

Solution *RandomizeTourHamiltonPath(const Matrix *distance_matrix, CustomVector<double> *input_tour) {
    CustomVector<double> *out_tour = new CustomVector<double>(*input_tour);
    out_tour->Shuffle(1, out_tour->GetSize() - 1);
    Solution *res = new Solution();
    res->tour = out_tour;
    res->distance = TourDistanceCalc(distance_matrix, out_tour);
    return res;
}

bool CheckDubl(CustomVector<double> *input_tour, int end = 0) {
    for (int i = 0; i < input_tour->GetSize()-1-end; ++i) {
        for (int j = i+1; j < input_tour->GetSize()-end; ++j) {
            if ((*input_tour)[i] == (*input_tour)[j]) {
                std::cout << i << " " << j << " " << (*input_tour)[i] << std::endl;
                return true;
            }
        }
    }
    return false;
}


void LocalSearchHamiltonPath(const Matrix *distance_matrix, Solution *sol) {

    unsigned int sol_len = sol->tour->GetSize();

    Solution best_route = Solution(*sol);
    Solution seed = Solution(*sol);

    for (unsigned int i = 1; i < sol_len - 2; ++i) {
        for (unsigned int j = i + 1; j < sol_len - 1; ++j) {
            best_route.tour->Reverse(i, j);
            best_route.distance = TourDistanceCalc(distance_matrix, best_route.tour);
            if (best_route.distance < sol->distance) {
                *sol = best_route;
            }
            best_route = seed;
        }
    }
}

Solution *CrossoverHamiltonPath(const Matrix *distance_matrix,
                    list<Solution *> reference_list,
                    double reverse_prob = 0.5,
                    double scramble_prob = 0.3) {
    std::random_device generator;
    std::uniform_int_distribution<unsigned int> int_distribution(0, reference_list.size() - 1);
    std::uniform_int_distribution<unsigned int> int_distribution_2(1, distance_matrix->GetRows() - 2);
    std::uniform_real_distribution<double> double_distribution(0, 1);
    unsigned int ix = int_distribution(generator);
    unsigned int iy = int_distribution(generator);
    while (ix == iy)
        iy = int_distribution(generator);

    CustomVector<double> *parent1 = new CustomVector<double>(distance_matrix->GetRows(), *GetElByIndex(reference_list, ix)->tour);
    CustomVector<double> *parent2 = new CustomVector<double>(distance_matrix->GetRows(), *GetElByIndex(reference_list, iy)->tour);


    unsigned int i = int_distribution_2(generator);
    unsigned int j = int_distribution_2(generator);
    while (i == j)
        j = int_distribution_2(generator);
    if (i > j)
        std::swap(i, j);

    double chance = double_distribution(generator);
    if (chance < reverse_prob) { // Прокнул реверс
        parent1->Reverse(i, j);
    }


    // В parent_2 остаются только те элементы, которых нет во фрагменте parent_1
    unsigned int insert_index = 1;
    bool can_insert;
    for (unsigned int k = 1; k < parent2->GetSize(); ++k) {
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
    if (chance < scramble_prob) {
        if (insert_index-1 != 1) {
            parent2->Shuffle(1, insert_index-1);
        }
    }
    insert_index = 1;
    for (unsigned int k = 1; k < i; ++k) {
        (*parent1)[k] = (*parent2)[insert_index];
        insert_index++;
    }
    for (unsigned int k = j + 1; k < parent1->GetSize(); ++k) {
        (*parent1)[k] = (*parent2)[insert_index];
        insert_index++;
    }

    delete parent2;
    Solution *sol = new Solution();
    sol->tour = parent1;
    sol->distance = TourDistanceCalc(distance_matrix, parent1);
    return sol;
}

Solution *StartHamiltonPath(
        bool input_debug,
        const Matrix *dist_matrix,
        unsigned int _iterations = 150,
        unsigned int _reference_size = 10,
        double _reverse_prob = 0.5,
        double _scramble_prob = 0.3,
        int start_ = 0,
        int end_ = 0)
        {
    unsigned int iterations = _iterations;      /*кол-во итераций*/
    unsigned int reference_size = _reference_size;   /*размер реф сета*/
    double reverse_prob = _reverse_prob;          /*вероятность переворота последовательности*/
    double scramble_prob = _scramble_prob;         /*вероятность перемешивания*/
    bool debug = input_debug;           /* Дебаг консоль*/
    int start = start_;
    int end = end_;


    if (start > end) {
        std::swap(start, end);
    }

    if (start == end) {
        throw std::invalid_argument("End must be greater than start");
    }

    unsigned int city_num = dist_matrix->GetRows();
    CustomVector<double> *city_tour = CreateTourHamiltonPath(dist_matrix,start,end);

    // std::cout << *city_tour << std::endl;

    double dist = TourDistanceCalc(dist_matrix, city_tour);

    unsigned int iter = 0;
    Solution *best = new Solution(city_tour, dist);
    list<Solution *> reference_list;

    for (int i = 0; i < reference_size; ++i) {
        auto tmp = RandomizeTourHamiltonPath(dist_matrix, city_tour);
        reference_list.push_back(tmp);
    }

    while (iter < iterations) {
        if (debug && iter % 15 == 0) {
            cout << "Iteration: " << iter + 1 << " Best sol: " << best->distance << endl;
        }

        iter++;
        list<Solution *> candidate_list;
        for (unsigned int i = 0; i < reference_size; ++i) {
            auto tmp = CrossoverHamiltonPath(dist_matrix, reference_list, reverse_prob, scramble_prob);
            candidate_list.push_back(tmp);
        }
        for (unsigned int i = 0; i < reference_size; ++i) {
            Solution *new_sol = GetElByIndex(candidate_list, i);
            LocalSearchHamiltonPath(dist_matrix, new_sol);

        }
        for (int i = 0; i < candidate_list.size(); ++i) {
            auto obj = candidate_list.front();
            reference_list.push_back(obj);
            candidate_list.pop_front();
        }
        reference_list.sort(CompareSolutions);
        unsigned int cnt = reference_list.size() - reference_size;
        for (unsigned int i = 0; i < cnt; ++i) {
            delete reference_list.back();
            reference_list.pop_back();
        }
        for (int i = 0; i < reference_size; ++i) {
            auto sol = GetElByIndex(reference_list, i);
            if (sol->distance < best->distance) {
                *best = *sol;
            }
        }

        while (!candidate_list.empty()) {
            delete candidate_list.back();
            candidate_list.pop_back();
        }
    }

    while (!reference_list.empty()) {
        delete reference_list.back();
        reference_list.pop_back();
    }

    delete city_tour;
    return best;
}

