#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <chrono>
#include <random>
#include <algorithm>

namespace structures {
    // Класс-наследник std::vector
    template<typename T>
    class CustomVector : public std::vector<T> {
    public:
        using std::vector<T>::vector; // Наследуем конструкторы базового класса

        // Дополнительные конструкторы
        explicit CustomVector(unsigned int size = 0, const T *data = nullptr) : std::vector<T>(size) {
            if (data) {
                std::copy(data, data + size, this->begin());
            }
        }

        explicit CustomVector(unsigned int size, const CustomVector &v) : CustomVector(size, v.data()) {
        }

        CustomVector(const CustomVector &v) : std::vector<T>(v) {
        } // Копирование

        explicit CustomVector(std::vector<T> &&v) noexcept : std::vector<T>(std::move(v)) {
        }

        void Shuffle(unsigned int start, unsigned int end) {
            if (start >= this->size() || end > this->size() || start >= end) {
                throw std::out_of_range("Invalid start or end indices for Shuffle.");
            }

            std::random_device rd;
            std::mt19937 g(rd());
            std::shuffle(this->begin() + start, this->begin() + end, g);
        }

        void Reverse(unsigned int start, unsigned int end) {
            if (start >= this->size() || end > this->size()) {
                throw std::out_of_range("Invalid start or end indices for Reverse.");
            }

            if (start > end) {
                std::swap(start, end);
            }

            unsigned int range = end - start;
            for (unsigned int k = 0; k < static_cast<unsigned int>(std::ceil(range / 2.0)); ++k) {
                std::swap((*this)[start + k], (*this)[end - k - 1]);
            }
        }

        friend std::ostream& operator<<(std::ostream& in, const CustomVector& v) {
            for (int i = 0; i < v.GetSize(); i++) {
                in << v[i] << "\t";
            }
            return in;
        }

        unsigned int GetSize() const {
            return static_cast<unsigned int>(this->size());
        }

        [[nodiscard]] std::vector<int> ToIntVector(int offset = 0) const {
            std::vector<int> result;
            result.reserve(this->size());

            if constexpr (std::is_same_v<T, int>) {
                // Если T == int, просто копируем
                result.assign(this->begin(), this->end());
            } else if constexpr (std::is_convertible_v<T, int>) {
                // Если T можно преобразовать в int (например, float, double)
                std::transform(
                    this->begin(),
                    this->end(),
                    std::back_inserter(result),
                    [offset](const T& val) { return static_cast<int>(val+offset); }
                );
            } else {
                // Ошибка компиляции, если T нельзя преобразовать в int
                static_assert(
                    std::is_convertible_v<T, int>,
                    "Cannot convert element type to int"
                );
            }

            return result;
        }

    };
}
#endif //VECTOR_H