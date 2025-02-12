#pragma once
#include <vector>
#include <initializer_list>
#include "densematrices.hpp"  // Assumes this file defines DenseMatrix<T> and: using int64 = long long;

template<typename T>
class Matrix;

//---------------------------------------------------------------------
// SparseMatrix class declaration using CRS (Compressed Row Storage)
// and CCS (Compressed Column Storage) formats.
//---------------------------------------------------------------------
template <typename T>
class SparseMatrix {
public:
    // Compressed Row Storage (CRS)
    struct CompressedRowStorage {
        std::vector<T> values;
        std::vector<int64> col_ind;
        std::vector<int64> row_ptr;

        // Set or update the element at (i,j).
        void set(const T& value, const int64& i, const int64& j);
        // Get (by reference) the element at (i,j). If not found, returns a reference to a static zero.
        T& get(const int64& i, const int64& j);
        // Overloaded function call operators.
        T& operator()(const int64& i, const int64& j);
        const T& operator()(const int64& i, const int64& j) const;
    };

    // Compressed Column Storage (CCS)
    struct CompressedColumnStorage {
        std::vector<T> values;
        std::vector<int64> row_ind;
        std::vector<int64> col_ptr;

        // Set or update the element at (i,j) (here i is the column index).
        void set(const T& value, const int64& i, const int64& j);
        // Get (by reference) the element at (i,j). If not found, returns a reference to a static zero.
        T& get(const int64& i, const int64& j);
        // Overloaded function call operators.
        T& operator()(const int64& i, const int64& j);
        const T& operator()(const int64& i, const int64& j) const;
    };

    // Constructors
    SparseMatrix() = default;
    SparseMatrix(const std::vector<T>& c_matrix, const int64& i, const int64& j);
    SparseMatrix(const std::vector<std::vector<T>>& matrix);
    SparseMatrix(const DenseMatrix<T>& matrix);
    SparseMatrix(const std::initializer_list<std::initializer_list<T>>& matrix);

    // Assignment operators
    SparseMatrix<T>& operator=(const SparseMatrix& other) = default;
    SparseMatrix<T>& operator=(const std::vector<std::vector<T>>& matrix);

    // Build functions for CRS
    CompressedRowStorage build_crs(const std::vector<T>& c_matrix, const int64& i, const int64& j) const;
    CompressedRowStorage build_crs(const std::vector<std::vector<T>>& matrix) const;
    CompressedRowStorage build_crs(const DenseMatrix<T>& matrix) const;
    CompressedRowStorage build_crs(const std::initializer_list<std::initializer_list<T>>& matrix);

    // Build functions for CCS
    CompressedColumnStorage build_ccs(const std::vector<T>& c_matrix, const int64& i, const int64& j) const;
    CompressedColumnStorage build_ccs(const std::vector<std::vector<T>>& matrix) const;
    CompressedColumnStorage build_ccs(const DenseMatrix<T>& matrix) const;
    CompressedColumnStorage build_ccs(const std::initializer_list<std::initializer_list<T>>& matrix);

    // Convert the sparse matrix to a dense 2D vector.
    std::vector<std::vector<T>> to_dense() const;

    CompressedRowStorage get_crs() const { return crs; };
    CompressedRowStorage get_ccs() const { return ccs; };

    T& operator()(const int64& i, const int64& j);
    const T& operator()(const int64& i, const int64& j) const;

private:
    mutable CompressedRowStorage crs;
    mutable CompressedColumnStorage ccs;

    friend class Matrix<T>;
    template <typename U>
    friend Matrix<U> add_sparse(const Matrix<U>&, const Matrix<U>&);    
};

#include "../src/smatreal.hpp"
