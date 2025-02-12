#pragma once
#include "densematrices.hpp"
#include "sparsematrices.hpp"  // This file provides the minimal SparseMatrix stub.
#include <iostream>
#include <stdexcept>
#include <algorithm>    // For std::min

// --- Constructor from a contiguous 1D vector ---
template <typename T>
DenseMatrix<T>::DenseMatrix(const std::vector<T>& c_matrix, const int64& i, const int64& j) 
{
    rows_ = i;
    cols_ = j;
    if (rows_ * cols_ != c_matrix.size()) 
    {
        std::cerr << "Error: The provided matrix vector has incorrect size.\n";
        throw std::invalid_argument("Incorrect matrix vector size");
    }
    matrix_ = c_matrix;
}

// --- Constructor from a 2D vector ---
template <typename T>
DenseMatrix<T>::DenseMatrix(const std::vector<std::vector<T>>& matrix) 
{
    // Validate if the provided 2D vector is empty or it's row is empty
    if (matrix.empty() || matrix[0].empty()) 
    {
        rows_ = cols_ = 0;
        std::cerr << "Error: Matrix rows or columns are empty.\n";
        return;
    }

    // Verify that all rows are of equal length.
    if (matrix.size() > 1)
    {
        for (uint64 i = 1; i < matrix.size(); ++i)
        {
            if (matrix[0].size() != matrix[i].size()) 
            {
                std::cerr << "Invalid matrix: row " << i << " size is " << matrix[i].size() 
                          << " instead of " << matrix[0].size() << ".\n";
                throw std::invalid_argument("Invalid matrix: rows have inconsistent sizes.");
            }
        }
    }

    // Calculate the total number of the elements 
    int64 total_size = static_cast<int64>(matrix.size()) * matrix[0].size();
    uint64 index = 0;

    matrix_.reserve(total_size);

    // Copy each row into the internal 1D storage
    for (const auto& row : matrix) 
    {
        // Since all rows have the same size, iterate over each value in the row.
        for (const auto& value : row) 
        {
            matrix_[index++] = value;
        }
    }

    rows_ = matrix.size();
    cols_ = matrix[0].size();
}

// --- Get a row as a vector ---
template <typename T>
std::vector<T> DenseMatrix<T>::row(const int64& i) const 
{
    // Check if the row index is valid
    if (i < 0 || i >= static_cast<int64>(rows_))
        throw std::out_of_range("Row index out of range.");

    // std::vector<T> row_result;
    // row_result.reserve(cols_);

    // for (uint64 j = 0; j < cols_; ++j) 
    // {
    //     row_result.push_back(matrix_[i * cols_ + j]);
    // }

    // return row_result;

    // Calculate the starting iterator for the row.
    auto start_iter = matrix_.begin() + i * cols_;

    // Create and return a vector from the contiguous range [start, start + cols_).
    // std::next(start_iter, cols_) doing next things:
    // - start_iter - pointer to the first element of the vector
    // - std::next() moves iterator on next cols_ steps , so it points on the
    // first element of the next row
    // std::vector<T>(iter_start, iter_end) - range of the iterators, 
    // copying elements between them
    return std::vector<T>(start_iter, std::next(start_iter, cols_));
}

// --- Get a column as a vector ---
template <typename T>
std::vector<T> DenseMatrix<T>::col(const int64& j) const 
{
    // Check if the column index is valid
    if (j < 0 || j >= static_cast<int64>(cols_))
        throw std::out_of_range("Column index out of range.");

    std::vector<T> col_result;
    col_result.reserve(rows_);

    for (uint64 i = 0; i < rows_; ++i) 
    {
        col_result[i] = matrix_[i * cols_ + j];
    }

    return col_result;
}

// --- Return number of rows ---
template <typename T>
inline int64 DenseMatrix<T>::row_count() const 
{
    return static_cast<int64>(rows_);
}

// --- Return number of columns ---
template <typename T>
inline int64 DenseMatrix<T>::col_count() const 
{
    return static_cast<int64>(cols_);
}

// --- Return the internal 1D matrix vector ---
template <typename T>
inline std::vector<T> DenseMatrix<T>::data() const 
{
    return matrix_;
}

// --- Convert to 2D vector ---
template <typename T>
std::vector<std::vector<T>> DenseMatrix<T>::mat() const 
{
    std::vector<std::vector<T>> result(rows_, std::vector<T>(cols_, T()));

    for (uint64 i = 0; i < rows_; ++i)
    {
        for (uint64 j = 0; j < cols_; ++j) 
        {
            result[i][j] = matrix_[i * cols_ + j];
        }
    }

    return result;
}

// --- Constructor from a SparseMatrix ---
template <typename T>
DenseMatrix<T>::DenseMatrix(const SparseMatrix<T>& sparse_m) 
{
    // The number of rows is inferred from the size of the CRS row pointer array
    rows_ = sparse_m.get_crs().row_ptr.size() - 1;

    // Convert sparse matrix to a dense 2D array
    std::vector<std::vector<T>> dense_from_sparse = sparse_m.to_dense();

    if (!dense_from_sparse.empty() && !dense_from_sparse[0].empty())
         cols_ = dense_from_sparse[0].size();
    else
         cols_ = 0;
    
    // Resize the internal storage to hold the matrix, initializing with default-constructed T.
    matrix_.resize(rows_ * cols_, T());

    // For each row, iterate over its nonzero elements and assign the corresponding value 
    for (int64 i = 0; i < static_cast<int64>(rows_); ++i) 
    {
        int64 start_index = sparse_m.get_crs().row_ptr[i];
        int64 end_index = sparse_m.get_crs().row_ptr[i + 1];

        // Loop over nonzero elements in row i
        for (int64 j = start_index; j < end_index; ++j) 
        {
            // Column index for current nonzero element 
            int64 col_index = sparse_m.get_crs().col_ind[j];
            // Place the value for correct position in the 1D vector representation
            matrix_[i * cols_ + col_index] = sparse_m.get_crs().values[j];
        }
    }
}

// --- Element access operator (non-const) ---
template <typename T>
T& DenseMatrix<T>::operator()(const int64& i, const int64& j) 
{
    if(i < 0 || i >= static_cast<int64>(rows_) || j < 0 || j >= static_cast<int64>(cols_))
        throw std::out_of_range("Index out of range.");

    return matrix_[i * cols_ + j];
}

// --- Element access operator (const) ---
template <typename T>
const T& DenseMatrix<T>::operator()(const int64& i, const int64& j) const 
{
    if(i < 0 || i >= static_cast<int64>(rows_) || j < 0 || j >= static_cast<int64>(cols_))
        throw std::out_of_range("Index out of range.");

    return matrix_[i * cols_ + j];
}

// --- Copy constructor ---
template <typename T>
DenseMatrix<T>::DenseMatrix(const DenseMatrix<T>& other) 
{
    std::cout << "Copy constructor called" << std::endl;
    matrix_ = other.matrix_;
    rows_ = other.rows_;
    cols_ = other.cols_;
    std::cout << "Copy constructor successful\n";
}

// --- Assignment operator from a 2D vector ---
template <typename T>
DenseMatrix<T>& DenseMatrix<T>::operator=(const std::vector<std::vector<T>>& matrix2d) 
{
    *this = DenseMatrix(matrix2d);
    return *this;
}

// --- Resize the matrix ---
template <typename T>
void DenseMatrix<T>::resize(const int64& rs, const int64& cs, const T& def_val) 
{
    // Allocate a new 1D vector with new dimensions, initializing all with def_val
    std::vector<T> new_matrix(rs * cs, def_val);

    // Determine the overlapping region dimensions.
    int64 min_rows = std::min(static_cast<int64>(rows_), rs);
    int64 min_cols = std::min(static_cast<int64>(cols_), cs);

    // for (int64 i = 0; i < min_rows; ++i) 
    // {
    //     for (int64 j = 0; j < min_cols; ++j) 
    //     {
    //         new_matrix[i * cs + j] = matrix_[i * cols_ + j];
    //     }
    // }

    for (int64 i = 0; i < min_rows; ++i) 
    {
        auto src_begin = matrix_.begin() + i * cols_;
        auto src_end   = src_begin + min_cols;
        auto dst_begin = new_matrix.begin() + i * cs;

        std::copy(src_begin, src_end, dst_begin); 
    }

    matrix_ = std::move(new_matrix);
    rows_ = rs;
    cols_ = cs;
}
