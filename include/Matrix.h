#pragma once
#include <ostream>

namespace my_matrix {
	class Row {
		unsigned int cols;
		double* row_data;
	public:
		Row(unsigned int cols, double* data) :cols(cols), row_data(data) {}
		~Row();
		double& operator[](unsigned int cols_);
		const double& operator[](unsigned int cols_) const;
	};

	class Matrix {
		unsigned int rows_;
		unsigned int cols_;
		double* data_;
		unsigned int id_;
		static unsigned int next_id;

		bool Can_Sum(const Matrix& other) const;
		bool Can_Mul(const Matrix& other) const;

	public:

		Matrix();
		Matrix(const unsigned int rows, const unsigned int cols, const double* data = nullptr);
		Matrix(const unsigned int n, const double* data = nullptr);
		Matrix(const Matrix& m);
		Matrix(Matrix&& m) noexcept;
		~Matrix();

		unsigned int Get_ID() const;
		unsigned int GetRows() const;
		unsigned int GetCols() const;


		Matrix& operator = (const Matrix& m);
		Matrix& operator = (Matrix&& m) noexcept;

		Matrix& operator += (const Matrix& m);
		Matrix& operator -= (const Matrix& m);
		Matrix& operator *= (const Matrix& m);

		Matrix& operator *= (const double scalar);

		const Row operator[](unsigned int rows_) const;
		Row operator[](unsigned int rows_);

		friend std::ostream& operator<<(std::ostream& in, const Matrix& m);
	};
	Matrix operator+(const Matrix& m1, const Matrix& m2);
	Matrix operator-(const Matrix& m1, const Matrix& m2);
	Matrix operator*(const Matrix& m1, const Matrix& m2);

	Matrix operator*(const Matrix& m1, const double scalar);
}
