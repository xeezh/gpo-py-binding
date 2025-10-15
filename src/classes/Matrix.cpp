#include <iostream>
#include <string>
#include "../../include/Matrix.h"

using namespace my_matrix;
using namespace std;

unsigned int Matrix::next_id = 1;

// class row
Row:: ~Row() {
	cols = 0;
	row_data = nullptr;
}
double& Row::operator[](unsigned int cols_) {
	if (cols_ > cols)
		throw std::invalid_argument("Index out of range");
	return row_data[cols_];
}

const double& Row::operator[](unsigned int cols_) const {
	if (cols_ > cols)
		throw std::invalid_argument("Index out of range");
	return row_data[cols_];
}

Matrix::Matrix() : rows_(0), cols_(0), data_(nullptr), id_(next_id++) {}

Matrix::Matrix(const unsigned int rows, const unsigned int cols, const double* data) {
	rows_ = rows;
	cols_ = cols;
	id_ = next_id++;
	if (rows*cols == 0) {
		data_ = nullptr;
		return;
	}
	data_ = new double[rows * cols];
	if (data == nullptr) {
		for (int i = 0; i < rows * cols; i++)
			data_[i] = 0;
	}
	else
		for (int i = 0; i < rows * cols; i++) data_[i] = data[i];
}

Matrix::Matrix(const unsigned int n, const double* data) : Matrix(n, n, data) {}

Matrix::Matrix(const Matrix& m) : Matrix(m.rows_, m.cols_, m.data_) {}

Matrix::Matrix(Matrix&& m) noexcept {
	swap(rows_, m.rows_);
	swap(cols_, m.cols_);
	swap(data_, m.data_);
}

Matrix::~Matrix() {
	cols_ = 0;
	rows_ = 0;
	if (data_ != nullptr)
		delete[] data_;
}

unsigned int Matrix::Get_ID() const {
	return id_;
}

unsigned int Matrix::GetRows() const {
	return rows_;
}

unsigned int Matrix::GetCols() const {
	return cols_;
}

bool Matrix::Can_Sum(const Matrix& other) const {
	return rows_ == other.rows_ && cols_ == other.cols_;
}

bool Matrix::Can_Mul(const Matrix& other) const {
	return cols_ == other.rows_;
}

Matrix& Matrix::operator=(const Matrix& m) {
	if (this != &m) {
		if (rows_ * cols_ != m.rows_ * m.cols_) {
			delete[] data_;
			data_ = new double[m.rows_ * m.cols_];
		}
		rows_ = m.rows_;
		cols_ = m.cols_;
		copy(m.data_, m.data_ + rows_ * cols_, data_);
	}
	return *this;
}

Matrix& Matrix::operator = (Matrix&& m) noexcept {
	if (this != &m) {
		swap(rows_, m.rows_);
		swap(cols_, m.cols_);
		swap(data_, m.data_);
	}
	return *this;
}

Matrix& Matrix::operator += (const Matrix& m) {
	if (!Can_Sum(m)) {
		string text = "Can't be sum " + to_string(this->Get_ID()) + " " + to_string(m.Get_ID());
		throw invalid_argument(text);
	}
	for (int i = 0; i < rows_ * cols_; ++i) {
		data_[i] += m.data_[i];
	}
	return *this;
}

Matrix& Matrix::operator -= (const Matrix& m) {
	if (!Can_Sum(m)) {
		string text = "Can't be sub " + to_string(this->Get_ID()) + " " + to_string(m.Get_ID());
		throw invalid_argument(text);
	}
	for (int i = 0; i < rows_ * cols_; ++i) {
		data_[i] -= m.data_[i];
	}
	return *this;
}

Matrix& Matrix::operator *= (const Matrix& m) {
	if (!Can_Mul(m)) {
		string text = "Can't be mul " + to_string(this->Get_ID()) + " " + to_string(m.Get_ID());
		throw invalid_argument(text);
	}
	for (int i = 0; i < rows_; i++) {
		for (int j = 0; j < m.cols_; j++) {
			double sum = 0;
			for (int k = 0; k < cols_; k++) {
				sum += data_[cols_ * i + k] * m[k][j];
			}
			(*this)[i][j] = sum;
		}
	}
	return *this;
}

Matrix& Matrix::operator *= (const double scalar) {
	for (int i = 0; i < rows_ * cols_; i++)
		data_[i] *= scalar;
	return *this;
}


Row Matrix::operator[](unsigned int rows_) {
	if (rows_ < 0 || rows_ >= this->rows_) {
		string text = "Index out of range " + to_string(this->Get_ID());
		throw invalid_argument(text);
	}
	Row current_row(cols_, &data_[rows_ * cols_]);
	return current_row;
}

const Row Matrix::operator[](unsigned int rows_) const {
	if (rows_ < 0 || rows_ >= this->rows_) {
		string text = "Index out of range " + to_string(this->Get_ID()) ; 
		throw invalid_argument(text);
	}
	Row current_row(cols_, &data_[rows_ * cols_]);
	return current_row;
}

ostream& my_matrix::operator<<(ostream& in, const Matrix& m) {
	for (int i = 0; i < m.GetRows(); i++) {
		for (int j = 0; j < m.GetCols(); j++) {
			in << m[i][j] << "\t";
		}
		in << endl;
	}
	return in;
}

Matrix my_matrix::operator+(const Matrix& m1, const Matrix& m2) {
	Matrix temp(m1);
	temp += m2;
	return temp;
}

Matrix my_matrix::operator-(const Matrix& m1, const Matrix& m2) {
	Matrix temp(m1);
	temp -= m2;
	return temp;
}

Matrix my_matrix::operator*(const Matrix& m1, const Matrix& m2) {
	Matrix temp(m1);
	temp *= m2;
	return temp;
}

Matrix my_matrix::operator*(const Matrix& m1, const double scalar) {
	Matrix temp(m1);
	temp *= scalar;
	return temp;
}