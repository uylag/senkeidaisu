#pragma once
#include <initializer_list>
#include "densematrices.hpp"
#include "sparsematrices.hpp"

template <typename T>
class Matrix {
private:
    SparseMatrix<T> sparse_;
    DenseMatrix<T> dense_;
    
protected:

public:
    Matrix() = default;
    Matrix(const std::vector<T>& c_matrix, const long long& i, const long long& j);
    Matrix(const std::vector<std::vector<T>>& matrix);
    Matrix(const std::initializer_list<T>& c_matrix);
    Matrix(const std::initializer_list<std::initializer_list<T>>& matrix);

    Matrix<T> operator+(const Matrix<T>& other);
    Matrix<T> operator+(const std::vector<std::vector<T>>& matrix);

    T& operator()(const int64& i, const int64& j);
    const T& operator()(const int64& i, const int64& j) const;

    friend Matrix<T> add_sparse<>(const Matrix<T>&, const Matrix<T>&);
};

#include "../src/matreal.hpp"