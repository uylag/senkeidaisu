#include "../include/densematrices.hpp"
#include <iostream>
#include <exception>

template <typename T>
DenseMatrix<T>::DenseMatrix() {
    matrix = std::vector<T>(1, 0);
}

template <typename T>
DenseMatrix<T>::DenseMatrix(const std::vector<T>& c_matrix, const long long& i, const long long& j) {
    if (rows * columns != c_matrix.size()) {
        std::cerr << "The matrix is in the wrong format.\n";
        return;
    }

    matrix = c_matrix;
    rows = i;
    columns = j;
}

template <typename T>
DenseMatrix<T>::DenseMatrix(const std::vector<std::vector<T>>& matrix) {
    if (matrix.empty() || matrix[0].empty()) {
        rows = columns = 0;
        std::cerr << "Matrix rows either columns are empty.\n";
        return;
    }

    for (int i = 1; i < matrix.size(); ++i) {
        if (matrix[0].size() != matrix[i].size()) {
            std::cerr << "Invalid matrix, given matrix has not equivalent size " 
                      << matrix[i].size() 
                      << " instead of " 
                      << matrix[0].size() << std::endl;
            throw std::invalid_argument("Invalid matrix: rows have inconsistent sizes.");
        }
    }

    long long total_size = matrix.size() * matrix[0].size();
    this->matrix.reserve(total_size);

    for (int i = 0; i < matrix.size(); ++i) {
        this->matrix.insert(this->matrix.end(), 
                            std::make_move_iterator(matrix[i].begin()), 
                            std::make_move_iterator(matrix[i].end()));
    }

    rows = matrix.size();
    columns = matrix[0].size();
}

template <typename T>
std::vector<T> DenseMatrix<T>::row(long long i) const {
    std::vector<T> row_result;
    for (int jdx = 0; jdx < columns; ++jdx) {
        row_result.push_back(matrix[i * columns + jdx]);
    }
    return row_result;
}

template <typename T>
std::vector<T> DenseMatrix<T>::col(long long j) const {
    std::vector<T> col_result;
    for (int idx = 0; idx < rows; ++idx) {
        col_result.push_back(matrix[idx * columns + j]);
    }
    return col_result;
}

template <typename T>
inline long long DenseMatrix<T>::row_count() const {
    return rows;
}

template <typename T>
inline long long DenseMatrix<T>::col_count() const {
    return columns;
}

template <typename T>
inline std::vector<T> DenseMatrix<T>::mat() const {
    return matrix;
}

