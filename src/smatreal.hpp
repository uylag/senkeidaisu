#pragma once
#include "../include/sparsematrices.hpp"
#include "../include/densematrices.hpp"
#include <iostream>
#include <algorithm>
#include <vector>
#include <stdexcept>

//---------------------------------------------------------------------
// Build CRS from a 2D dense matrix.
// For each row, stores nonzero elements and their column indices.
//---------------------------------------------------------------------
template <typename T>
typename SparseMatrix<T>::CompressedRowStorage SparseMatrix<T>::build_crs(const std::vector<std::vector<T>>& matrix) const 
{
    typename SparseMatrix<T>::CompressedRowStorage crs_result;
    crs_result.row_ptr.push_back(0);
    for (size_t i = 0; i < matrix.size(); ++i) {
        for (size_t j = 0; j < matrix[i].size(); ++j) {
            if (matrix[i][j] != 0) {
                crs_result.values.push_back(matrix[i][j]);
                crs_result.col_ind.push_back(static_cast<int64>(j));
            }
        }
        crs_result.row_ptr.push_back(static_cast<int64>(crs_result.values.size()));
    }
    crs = crs_result;
    return crs_result;
}

//---------------------------------------------------------------------
// Build CRS from a 1D contiguous matrix.
// 'i' is the number of rows and 'j' the number of columns.
//---------------------------------------------------------------------
template <typename T>
typename SparseMatrix<T>::CompressedRowStorage SparseMatrix<T>::build_crs(const std::vector<T>& c_matrix, const int64& i, const int64& j) const
{
    typename SparseMatrix<T>::CompressedRowStorage crs_result;
    if (i == 0 || j == 0)
        return crs_result;
    crs_result.row_ptr.push_back(0);
    for (int64 row = 0; row < i; ++row) {
        for (int64 col = 0; col < j; ++col) {
            int64 index = row * j + col;
            if (index >= static_cast<int64>(c_matrix.size())) {
                std::cerr << "Matrix given in wrong form.\n";
                return crs_result;
            }
            if (c_matrix[index] != 0) {
                crs_result.values.push_back(c_matrix[index]);
                crs_result.col_ind.push_back(col);
            }
        }
        crs_result.row_ptr.push_back(static_cast<int64>(crs_result.values.size()));
    }
    crs = crs_result;
    return crs_result;
}

//---------------------------------------------------------------------
// Build CRS from a DenseMatrix by extracting its 1D data vector and dimensions.
//---------------------------------------------------------------------
template <typename T>
typename SparseMatrix<T>::CompressedRowStorage SparseMatrix<T>::build_crs(const DenseMatrix<T>& matrix) const{
    return build_crs(matrix.mat());
}

//---------------------------------------------------------------------
// Build CCS from a 2D dense matrix.
// Outer loop is over columns; for each column, stores nonzero elements and their row indices.
//---------------------------------------------------------------------
template <typename T>
typename SparseMatrix<T>::CompressedColumnStorage SparseMatrix<T>::build_ccs(const std::vector<std::vector<T>>& matrix) const {
    typename SparseMatrix<T>::CompressedColumnStorage ccs_result;
    if (matrix.empty() || matrix[0].empty())
        return ccs_result;
    
    int nrows = static_cast<int>(matrix.size());
    int ncols = static_cast<int>(matrix[0].size());
    
    ccs_result.col_ptr.push_back(0);
    for (int col = 0; col < ncols; ++col) {
        for (int row = 0; row < nrows; ++row) {
            if (matrix[row][col] != 0) {
                ccs_result.values.push_back(matrix[row][col]);
                ccs_result.row_ind.push_back(row);
            }
        }
        ccs_result.col_ptr.push_back(static_cast<int64>(ccs_result.values.size()));
    }
    ccs = ccs_result;
    return ccs_result;
}

//---------------------------------------------------------------------
// Build CCS from a 1D contiguous matrix.
// The matrix is interpreted in row‑major order; we iterate by columns to build CCS.
//---------------------------------------------------------------------
template <typename T>
typename SparseMatrix<T>::CompressedColumnStorage SparseMatrix<T>::build_ccs(const std::vector<T>& c_matrix, const int64& i, const int64& j) const
{
    typename SparseMatrix<T>::CompressedColumnStorage ccs_result;
    if (i == 0 || j == 0)
        return ccs_result;
    int64 nrows = i;
    int64 ncols = j;
    ccs_result.col_ptr.push_back(0);
    for (int64 col = 0; col < ncols; ++col) {
        for (int64 row = 0; row < nrows; ++row) {
            int64 index = row * j + col;
            if (index >= static_cast<int64>(c_matrix.size())) {
                std::cerr << "Matrix given in wrong form.\n";
                return ccs_result;
            }
            if (c_matrix[index] != 0) {
                ccs_result.values.push_back(c_matrix[index]);
                ccs_result.row_ind.push_back(row);
            }
        }
        ccs_result.col_ptr.push_back(static_cast<int64>(ccs_result.values.size()));
    }
    ccs = ccs_result;
    return ccs_result;
}

//---------------------------------------------------------------------
// Build CCS from a DenseMatrix.
//---------------------------------------------------------------------
template <typename T>
typename SparseMatrix<T>::CompressedColumnStorage SparseMatrix<T>::build_ccs(const DenseMatrix<T>& matrix) const 
{
    return build_ccs(matrix.mat(), matrix.row_count(), matrix.col_count());
}

//---------------------------------------------------------------------
// Get element from CompressedRowStorage by row and column indices.
// If not found (i.e. zero), returns a reference to a static zero.
// (Modifying the returned static zero is not recommended.)
//---------------------------------------------------------------------
template <typename T>
T& SparseMatrix<T>::CompressedRowStorage::get(const int64& i, const int64& j) {
    static T zero = T();
    if (i < 0 || i >= static_cast<int64>(row_ptr.size() - 1))
        return zero;
    for (int64 idx = row_ptr[i]; idx < row_ptr[i + 1]; ++idx) {
        if (col_ind[idx] == j)
            return values[idx];
    }
    return zero;
}

template <typename T>
const T& SparseMatrix<T>::CompressedRowStorage::operator()(const int64& i, const int64& j) const {
    return const_cast<CompressedRowStorage*>(this)->get(i, j);
}

template <typename T>
T& SparseMatrix<T>::CompressedRowStorage::operator()(const int64& i, const int64& j) {
    return get(i, j);
}

//---------------------------------------------------------------------
// Get element from CompressedColumnStorage by column and row indices.
// Here, i is treated as the column index and j as the row index.
// If not found, returns a reference to a static zero.
//---------------------------------------------------------------------
template <typename T>
T& SparseMatrix<T>::CompressedColumnStorage::get(const int64& i, const int64& j) {
    static T zero = T();
    if (i < 0 || i >= static_cast<int64>(col_ptr.size() - 1))
        return zero;
    for (int64 idx = col_ptr[i]; idx < col_ptr[i + 1]; ++idx) {
        if (row_ind[idx] == j)
            return values[idx];
    }
    return zero;
}

template <typename T>
const T& SparseMatrix<T>::CompressedColumnStorage::operator()(const int64& i, const int64& j) const {
    return const_cast<CompressedColumnStorage*>(this)->get(i, j);
}

template <typename T>
T& SparseMatrix<T>::CompressedColumnStorage::operator()(const int64& i, const int64& j) {
    return get(i, j);
}

template <typename T>
const T& SparseMatrix<T>::operator()(const int64& i, const int64& j) const {
    return const_cast<CompressedColumnStorage*>(this)->crs.get(i, j);
}

template <typename T>
T& SparseMatrix<T>::operator()(const int64& i, const int64& j) {
    return crs.get(i, j);
}

//---------------------------------------------------------------------
// Set element in CompressedColumnStorage at column i and row j.
// If the element exists, update its value; otherwise, insert it and adjust col_ptr accordingly.
//---------------------------------------------------------------------
template <typename T>
void SparseMatrix<T>::CompressedColumnStorage::set(const T& value, const int64& i, const int64& j) {
    if (i < 0 || i >= static_cast<int64>(col_ptr.size() - 1))
        return;
    for (int64 idx = col_ptr[i]; idx < col_ptr[i + 1]; ++idx) {
        if (row_ind[idx] == j) {
            values[idx] = value;
            return;
        }
    }
    // Element not found; insert at the end of column i's segment.
    values.insert(values.begin() + col_ptr[i + 1], value);
    row_ind.insert(row_ind.begin() + col_ptr[i + 1], j);
    // Increment subsequent col_ptr entries.
    for (size_t k = i + 1; k < col_ptr.size(); ++k)
        col_ptr[k]++;
}

//---------------------------------------------------------------------
// Set element in CompressedRowStorage at row i and column j.
// If the element exists, update its value; otherwise, insert it and adjust row_ptr accordingly.
//---------------------------------------------------------------------
template <typename T>
void SparseMatrix<T>::CompressedRowStorage::set(const T& value, const int64& i, const int64& j) {
    if (i < 0 || i >= static_cast<int64>(row_ptr.size() - 1))
        return;
    for (int64 idx = row_ptr[i]; idx < row_ptr[i + 1]; ++idx) {
        if (col_ind[idx] == j) {
            values[idx] = value;
            return;
        }
    }
    // Element not found; insert at the end of row i's segment.
    values.insert(values.begin() + row_ptr[i + 1], value);
    col_ind.insert(col_ind.begin() + row_ptr[i + 1], j);
    // Increment subsequent row_ptr entries.
    for (size_t k = i + 1; k < row_ptr.size(); ++k)
        row_ptr[k]++;
}

//---------------------------------------------------------------------
// SparseMatrix constructor from a 1D contiguous vector.
//---------------------------------------------------------------------
template <typename T>
SparseMatrix<T>::SparseMatrix(const std::vector<T>& c_matrix, const int64& i, const int64& j) {
    build_crs(c_matrix, i, j);
    build_ccs(c_matrix, i, j);
}

//---------------------------------------------------------------------
// SparseMatrix constructor from a 2D vector.
//---------------------------------------------------------------------
template <typename T>
SparseMatrix<T>::SparseMatrix(const std::vector<std::vector<T>>& matrix) {
    build_crs(matrix);
    build_ccs(matrix);
}

//---------------------------------------------------------------------
// Convert the sparse matrix to a dense 2D vector.
// Determines the number of rows from CRS and the maximum column index from crs.col_ind.
//---------------------------------------------------------------------
template <typename T>
std::vector<std::vector<T>> SparseMatrix<T>::to_dense() const {
    int64 nrows = static_cast<int64>(crs.row_ptr.size()) - 1;
    int64 max_columns = 0;
    for (const auto& col : crs.col_ind) {
        max_columns = std::max(max_columns, col + 1);
    }
    std::vector<std::vector<T>> dense_matrix(nrows, std::vector<T>(max_columns, T()));
    for (int64 i = 0; i < nrows; ++i) {
        for (int64 idx = crs.row_ptr[i]; idx < crs.row_ptr[i + 1]; ++idx) {
            int64 col = crs.col_ind[idx];
            dense_matrix[i][col] = crs.values[idx];
        }
    }
    return dense_matrix;
}

//---------------------------------------------------------------------
// SparseMatrix constructor from a DenseMatrix.
//---------------------------------------------------------------------
template <typename T>
SparseMatrix<T>::SparseMatrix(const DenseMatrix<T>& matrix) {
    build_crs(matrix);
    build_ccs(matrix);
}

//---------------------------------------------------------------------
// Assignment operator from a 2D vector.
// Constructs a temporary SparseMatrix and moves it into *this.
//---------------------------------------------------------------------
template <typename T>
SparseMatrix<T>& SparseMatrix<T>::operator=(const std::vector<std::vector<T>>& matrix) {
    *this = SparseMatrix(matrix);
    return *this;
}

template <typename T>
typename SparseMatrix<T>::CompressedRowStorage
SparseMatrix<T>::build_crs(const std::initializer_list<std::initializer_list<T>>& matrix) {
    int64 rows, cols;
    if (matrix.size() == 1) {
        const auto& first = *matrix.begin();
        rows = std::distance(first.begin(), first.end());
        cols = 1;
    } else {
        rows = matrix.size();
        const auto& first = *matrix.begin();
        cols = std::distance(first.begin(), first.end());
    }
    CompressedRowStorage crs_result;
    crs_result.row_ptr.push_back(0);
    if (matrix.size() == 1) {
        int64 rowIndex = 0;
        for (auto it = matrix.begin()->begin(); it != matrix.begin()->end(); ++it) {
            if (*it != 0) {
                crs_result.values.push_back(*it);
                crs_result.col_ind.push_back(0);
            }
            crs_result.row_ptr.push_back(static_cast<int64>(crs_result.values.size()));
            rowIndex++;
        }
    } else {
        for (const auto& rowList : matrix) {
            int64 colIndex = 0;
            for (const auto& elem : rowList) {
                if (elem != 0) {
                    crs_result.values.push_back(elem);
                    crs_result.col_ind.push_back(colIndex);
                }
                colIndex++;
            }
            crs_result.row_ptr.push_back(static_cast<int64>(crs_result.values.size()));
        }
    }
    crs = crs_result;
    return crs_result;
}

template <typename T>
typename SparseMatrix<T>::CompressedColumnStorage
SparseMatrix<T>::build_ccs(const std::initializer_list<std::initializer_list<T>>& matrix) {
    int64 rows, cols;
    if (matrix.size() == 1) {
        const auto& first = *matrix.begin();
        rows = std::distance(first.begin(), first.end());
        cols = 1;
    } else {
        rows = matrix.size();
        const auto& first = *matrix.begin();
        cols = std::distance(first.begin(), first.end());
    }
    CompressedColumnStorage ccs_result;
    ccs_result.col_ptr.push_back(0);
    if (matrix.size() == 1) {
        int64 rowIndex = 0;
        for (auto it = matrix.begin()->begin(); it != matrix.begin()->end(); ++it) {
            if (*it != 0) {
                ccs_result.values.push_back(*it);
                ccs_result.row_ind.push_back(rowIndex);
            }
            rowIndex++;
        }
        ccs_result.col_ptr.push_back(static_cast<int64>(ccs_result.values.size()));
    } else {
        std::vector<std::initializer_list<T>> outer(matrix);
        int64 nrows = outer.size();
        int64 ncols = std::distance(outer[0].begin(), outer[0].end());
        for (int64 col = 0; col < ncols; ++col) {
            for (int64 row = 0; row < nrows; ++row) {
                std::vector<T> vec(outer[row].begin(), outer[row].end());
                if (vec[col] != 0) {
                    ccs_result.values.push_back(vec[col]);
                    ccs_result.row_ind.push_back(row);
                }
            }
            ccs_result.col_ptr.push_back(static_cast<int64>(ccs_result.values.size()));
        }
    }
    ccs = ccs_result;
    return ccs_result;
}

template <typename T>
SparseMatrix<T>::SparseMatrix(const std::initializer_list<std::initializer_list<T>>& matrix) {
    build_crs(matrix);
    build_ccs(matrix);
}
