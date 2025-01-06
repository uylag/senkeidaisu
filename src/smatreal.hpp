#pragma once
#include "../include/sparsematrices.hpp"
#include "../include/densematrices.hpp"
#include <iostream>
#include <algorithm>

template <typename T>
typename SparseMatrix<T>::CompressedRowStorage SparseMatrix<T>::build_crs(const std::vector<std::vector<T>>& matrix) {
    SparseMatrix<T>::CompressedRowStorage crs_result;

    crs_result.row_ptr.push_back(0);
    for (int i = 0; i < matrix.size(); ++i) {
        for (int j = 0; j < matrix[0].size(); ++j) {
            if (matrix[i][j] != 0) {
                crs_result.values.push_back(matrix[i][j]);
                crs_result.col_ind.push_back(j);
            }
        }
        crs_result.row_ptr.push_back(crs_result.values.size());
    }

    crs = crs_result;
    return crs_result;
}

template <typename T>
typename SparseMatrix<T>::CompressedRowStorage SparseMatrix<T>::build_crs(const std::vector<T>& c_matrix, const long long& i, const long long& j) {
    SparseMatrix<T>::CompressedRowStorage crs_result;

    if (i == 0 || j == 0)
        return crs_result;

    crs_result.row_ptr.push_back(0);
    for (int idx = 0; idx < i; ++idx) {
        for (int jdx = 0; jdx < j; ++jdx) {
            if (idx * j + jdx >= c_matrix.size()) {
                std::cerr << "Matrix given in wrong form.\n";
                return crs_result;
            }
            if (c_matrix[idx * j + jdx] != 0) {
                crs_result.values.push_back(c_matrix[idx * j + jdx]);
                crs_result.col_ind.push_back(jdx);
            }
        }
        crs_result.row_ptr.push_back(crs_result.values.size());
    }

    crs = crs_result;
    return crs_result;
}

template <typename T>
typename SparseMatrix<T>::CompressedRowStorage SparseMatrix<T>::build_crs(const DenseMatrix<T>& matrix) {
    SparseMatrix<T>::CompressedRowStorage crs_result;

    crs_result = build_crs(matrix.mat(), matrix.row_count(), matrix.col_count());
    return crs_result;
}

template <typename T>
typename SparseMatrix<T>::CompressedColumnStorage SparseMatrix<T>::build_ccs(const std::vector<std::vector<T>>& matrix) {
    SparseMatrix<T>::CompressedColumnStorage ccs_result;

    ccs_result.col_ptr.push_back(0);
    for (int i = 0; i < matrix.size(); ++i) {
        for (int j = 0; j < matrix[0].size(); ++j) {
            if (matrix[i][j] != 0) {
                ccs_result.values.push_back(matrix[i][j]);
                ccs_result.row_ind.push_back(j);
            }
        }
        ccs_result.col_ptr.push_back(ccs_result.values.size());
    }

    ccs = ccs_result;
    return ccs_result;
}

template <typename T>
typename SparseMatrix<T>::CompressedColumnStorage SparseMatrix<T>::build_ccs(const std::vector<T>& c_matrix, const long long& i, const long long& j) {
    SparseMatrix<T>::CompressedColumnStorage ccs_result;

    if (i == 0 || j == 0)
        return ccs_result;

    ccs_result.col_ptr.push_back(0);
    for (int idx = 0; idx < i; ++idx) {
        for (int jdx = 0; jdx < j; ++jdx) {
            if (idx * j + jdx >= c_matrix.size()) {
                std::cerr << "Matrix given in wrong form.\n";
                return ccs_result;
            }
            if (c_matrix[idx * j + jdx] != 0) {
                ccs_result.values.push_back(c_matrix[idx * j + jdx]);
                ccs_result.row_ind.push_back(jdx);
            }
        }
        ccs_result.col_ptr.push_back(ccs_result.values.size());
    }

    ccs = ccs_result;
    return ccs_result;
}

template <typename T>
typename SparseMatrix<T>::CompressedColumnStorage SparseMatrix<T>::build_ccs(const DenseMatrix<T>& matrix) {
    SparseMatrix<T>::CompressedColumnStorage ccs_result;

    ccs_result = build_ccs(matrix.mat(), matrix.row_count(), matrix.col_count());
    return ccs_result;
}

template <typename T>
T& SparseMatrix<T>::CompressedRowStorage::get(const long long& i, const long long& j) {
    static T zero = 0;

    if (i + 1 >= row_ptr.size()) {
        return zero;
    }

     for (int jdx = row_ptr[i]; jdx < row_ptr[i + 1]; ++jdx) {
        if (col_ind[jdx] == j) {
            return values[jdx];
        }
     }

     return zero;
}

template <typename T>
T& SparseMatrix<T>::CompressedColumnStorage::get(const long long& i, const long long& j) {
    static T zero = 0;

    if (i + 1 >= col_ptr.size()) {
        return zero;
    }

     for (int jdx = col_ptr[i]; jdx < col_ptr[i + 1]; ++jdx) {
        if (row_ind[jdx] == j) {
            return values[jdx];
        }
     }

     return zero;
}

template <typename T>
void SparseMatrix<T>::CompressedColumnStorage::set(const T& value, const long long& i, const long long& j) {
    if (i + 1 >= col_ptr.size()) {
        return;
    }

     for (int jdx = col_ptr[i]; jdx < col_ptr[i + 1]; ++jdx) {
        if (row_ind[jdx] == j) {
            values[jdx] = value;
            return;
        }
     }

    values.insert(values.begin() + col_ptr[i + 1], value);
    row_ind.insert(row_ind.begin() + col_ptr[i + 1], j);

    for (int k = i + 1; k < col_ptr.size(); ++k) {
        col_ptr[k]++;
    }
}

template <typename T>
void SparseMatrix<T>::CompressedRowStorage::set(const T& value, const long long& i, const long long& j) {
    if (i + 1 >= row_ptr.size()) {
        return;
    }

     for (int jdx = row_ptr[i]; jdx < row_ptr[i + 1]; ++jdx) {
        if (col_ind[jdx] == j) {
            values[jdx] = value;
            return;
        }
     }

    values.insert(values.begin() + row_ptr[i + 1], value);
    col_ind.insert(col_ind.begin() + row_ptr[i + 1], j);

    for (int k = i + 1; k < row_ptr.size(); ++k) {
        row_ptr[k]++;
    }
}

template <typename T>
const T& SparseMatrix<T>::CompressedRowStorage::operator()(const long long& i, const long long& j) const {
    return get(i, j);
}

template <typename T>
T& SparseMatrix<T>::CompressedRowStorage::operator()(const long long& i, const long long& j) {
    return get(i, j);
}

template <typename T>
const T& SparseMatrix<T>::CompressedColumnStorage::operator()(const long long& i, const long long& j) const {
    return get(i, j);
}

template <typename T>
T& SparseMatrix<T>::CompressedColumnStorage::operator()(const long long& i, const long long& j) {
    return get(i, j);
}

template <typename T>
SparseMatrix<T>::SparseMatrix(const std::vector<T>& c_matrix, const long long& i, const long long& j) {
    build_crs(c_matrix, i, j);
    build_ccs(c_matrix, i, j);
}

template <typename T>
SparseMatrix<T>::SparseMatrix(const std::vector<std::vector<T>>& matrix) {
    build_crs(matrix);
    build_ccs(matrix);
}