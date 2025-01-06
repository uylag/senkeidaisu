#pragma once
#include "../include/matrix.hpp"

template <typename T>
Matrix<T>::Matrix(const std::vector<T>& c_matrix, const long long& i, const long long& j) 
    : sparse(c_matrix, i, j), dense(c_matrix, i, j) {}

template <typename T>
Matrix<T>::Matrix(const std::vector<std::vector<T>> matrix)
     : sparse(matrix), dense(matrix) {}