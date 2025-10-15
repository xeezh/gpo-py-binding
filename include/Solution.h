#ifndef SOLUTION_H
#define SOLUTION_H

#include "Matrix.h"
#include "../src/classes/CustomVector.h"
#include <list>
#include <limits>
#include <algorithm>

using namespace my_matrix;

namespace structures {
    struct Solution {
        CustomVector<double>* tour;
        double distance;

        Solution();
        ~Solution();
        explicit Solution(const CustomVector<double>* other_tour, double other_distance);
        explicit Solution(const Solution& other);
        explicit Solution(Solution&& other);

        Solution& operator=(const Solution& other);
    };

    bool CompareSolutions(const Solution* s1, const Solution* s2);

    // Взятие элемента списка по индексу
    template<typename T>
    T& GetElByIndex(std::list<T>& lst, size_t index) {
        auto it = lst.begin();
        std::advance(it, index);
        return *it;
    }

    // Вычисление расстояния по туру
    template<typename T>
    double TourDistanceCalc(const Matrix* distance_matrix, const CustomVector<T>* city_tour) {
        double distance = 0;
        for (unsigned int i = 0; i < city_tour->GetSize() - 1; ++i) {
            unsigned int j = i + 1;
            distance += ((*distance_matrix)[(unsigned int)((*city_tour)[i]) - 1][(unsigned int)((*city_tour)[j]) - 1]);
        }
        return distance;
    }

    Solution* RandomizeTour(const Matrix* distance_matrix, CustomVector<double>* input_tour);
}

#endif