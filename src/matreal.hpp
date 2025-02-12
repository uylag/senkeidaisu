#pragma once
#include "../include/matrix.hpp"
#include <algorithm>
#include <initializer_list>
#include <iostream>

template <typename T>
Matrix<T>::Matrix(const std::vector<T>& c_matrix, const long long& i, const long long& j) 
    : sparse_(c_matrix, i, j), dense_(c_matrix, i, j) {}

template <typename T>
Matrix<T>::Matrix(const std::vector<std::vector<T>>& matrix)
    : sparse_(matrix), dense_(matrix) {}

template <typename T>
Matrix<T>::Matrix(const std::initializer_list<std::initializer_list<T>>& matrix)
    : sparse_(matrix), dense_(matrix) {};

template <typename T>
Matrix<T> add_sparse(const Matrix<T>& A, const Matrix<T>& B)
{
    // We need to ensure that both A and B have sparse representation
    // Get CRS representation from each
    auto A_crs = A.sparse_.build_crs(A.dense_);
    auto B_crs = B.sparse_.build_crs(B.dense_);

    // Number of rows is determined by the CRS row_ptr
    int64 nrows = static_cast<int64>(A_crs.row_ptr.size()) - 1;

    // We will build a new CRS for the result
    typename SparseMatrix<T>::CompressedRowStorage C_crs;
    C_crs.row_ptr.push_back(0);

    for (int64 i = 0; i < nrows; ++i)
    {
        int64 a_start = A_crs.row_ptr[i], a_end = A_crs.row_ptr[i + 1];
        int64 b_start = B_crs.row_ptr[i], b_end = B_crs.row_ptr[i + 1];
        int64 a_idx = a_start, b_idx = b_start;

        // Temporary vectors for row data
        std::vector<T> row_values;
        std::vector<int64> row_cols;

        // Merge the two lists until both lists are fully processed
        // - `a_idx < a_end` -> still elements left in A
        // - `b_idx < b_end` -> still elements left in B
        // - `||` instead of `&&` ensures we don't stop if one list is empty 
        while (a_idx < a_end || b_idx < b_end)
        {
            // Case 1: A's column index is smaller -> take element from A
            // - `a_idx < a_end` -> A still has elements
            // - `b_idx >= b_end` -> B is fully processed (take A)
            // - OR `A_crs.col_ind[a_idx] < B_crs.col_ind[b_idx]` → A's column comes first
            //
            //  - `A_crs.col_ind[a_idx]` -> Column index of current element in A
            //  - `B_crs.col_ind[b_idx]` -> Column index of current element in B
            //  - If `A_crs.col_ind[a_idx] < B_crs.col_ind[b_idx]`, then:
            //  A's element appears earlier in the row than B's element
            //  We must insert A's element first to maintain correct column order
            if (a_idx < a_end && (b_idx >= b_end || A_crs.col_ind[a_idx] < B_crs.col_ind[b_idx]))
            {
                row_values.push_back(A_crs.values[a_idx]);
                row_cols.push_back(A_crs.col_ind[a_idx]);
                a_idx++;
            }
            else if (b_idx < b_end && (a_idx >= a_end || B_crs.col_ind[b_idx] < A_crs.col_ind[a_idx]))
            {
                row_values.push_back(B_crs.values[b_idx]);
                row_cols.push_back(B_crs.col_ind[b_idx]);
                b_idx++;
            }
            else 
            {
                // Both matrices have a nonzero at this column.
                int64 col = A_crs.col_ind[a_idx]; // same as B_crs.col_ind[b_idx]
                T sum = A_crs.values[a_idx] + B_crs.values[b_idx];
                // Only store if the sum is nonzero.
                if (sum != 0)
                {
                    row_values.push_back(sum);
                    row_cols.push_back(col);
                }
                a_idx++;
                b_idx++;
            }
        }

        for (auto& val : row_values)
        {
            C_crs.values.push_back(val);
        }
        for (auto& col : row_cols)
        {
            C_crs.col_ind.push_back(col);
        }
        // New row pointer value
        C_crs.row_ptr.push_back(static_cast<int64>(C_crs.values.size()));
    }

    SparseMatrix<T> C;
    C.crs = C_crs;
    
    Matrix<T> result;
    result.sparse_ = C;
    result.dense_ = DenseMatrix<T>(C);

    return result;
}

template <typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& other)
{
    return add_sparse(*this, other);
}

template <typename T>
Matrix<T> Matrix<T>::operator+(const std::vector<std::vector<T>>& matrix)
{
    Matrix<T> other(matrix);
    return (*this) + other;
}

template <typename T>
T& Matrix<T>::operator()(const int64& i, const int64& j) 
{
    return sparse_(i, j);
};

template <typename T>
const T& Matrix<T>::operator()(const int64& i, const int64& j) const 
{
    return sparse_(i, j);
};