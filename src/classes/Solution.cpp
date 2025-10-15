#include "../../include/Solution.h"

namespace structures {
    // Конструктор по умолчанию
    Solution::Solution() {
        tour = nullptr;
        distance = std::numeric_limits<double>::max();
    }

    // Деструктор
    Solution::~Solution() {
        if (tour != nullptr)
            delete tour;
        distance = 0;
    }

    // Конструктор от других данных
    Solution::Solution(const CustomVector<double>* other_tour, double other_distance) {
        tour = new CustomVector<double>(*other_tour);
        distance = other_distance;
    }

    // Конструктор копирования
    Solution::Solution(const Solution& other) {
        tour = new CustomVector<double>(*other.tour);
        distance = other.distance;
    }

    // Конструктор перемещения
    Solution::Solution(Solution&& other) {
        tour = new CustomVector<double>();
        std::swap(tour, other.tour);
        distance = other.distance;
    }

    // Оператор присваивания
    Solution& Solution::operator=(const Solution& other) {
        if (tour == nullptr || tour->GetSize() != other.tour->GetSize()) {
            delete tour;
            tour = new CustomVector<double>(other.tour->GetSize());
        }
        for (int i = 0; i < tour->GetSize(); ++i) {
            (*tour)[i] = (*other.tour)[i];
        }
        distance = other.distance;
        return *this;
    }

    // Сравнение решений
    bool CompareSolutions(const Solution* s1, const Solution* s2) {
        return s1->distance < s2->distance;
    }

    // Функция генерации случайного маршрута
    Solution* RandomizeTour(const Matrix* distance_matrix, CustomVector<double>* input_tour) {
        CustomVector<double>* out_tour = new CustomVector<double>(*input_tour);
        out_tour->Shuffle(0, out_tour->GetSize() - 1);
        (*out_tour)[out_tour->GetSize() - 1] = (*out_tour)[0];
        Solution* res = new Solution();
        res->tour = out_tour;
        res->distance = TourDistanceCalc(distance_matrix, out_tour);
        return res;
    }
}