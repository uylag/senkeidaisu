#pragma once
#include <vector>
#include <initializer_list>

// Type aliases for clarity.
using int64 = long long;
using uint64 = unsigned long long;

// Forward declaration for SparseMatrix.
template <typename T> class SparseMatrix;

template <typename T>
class DenseMatrix {
private:
    uint64 rows_;
    uint64 cols_;
    std::vector<T> matrix_;

public:
    DenseMatrix() : rows_(0), cols_(0) {}

    // Constructs a DenseMatrix from a contiguous 1D vector.
    DenseMatrix(const std::vector<T>& c_matrix, const long long& i, const long long& j);

    // Constructs a DenseMatrix from a 2D vector.
    DenseMatrix(const std::vector<std::vector<T>>& matrix);

    // Constructs a DenseMatrix from a SparseMatrix.
    DenseMatrix(const SparseMatrix<T>& sparse_m);

    // Copy constructor.
    DenseMatrix(const DenseMatrix<T>& other);

    DenseMatrix(const std::initializer_list<std::initializer_list<T>>& matrix)
    {
        *this = DenseMatrix(SparseMatrix(matrix));
    }

    // Element access operators.
    T& operator()(const long long& i, const long long& j);
    const T& operator()(const long long& i, const long long& j) const;

    // Assignment operators.
    DenseMatrix<T>& operator=(const DenseMatrix& other) = default;
    DenseMatrix<T>& operator=(const std::vector<std::vector<T>>& matrix);

    // Get a specific row or column.
    std::vector<T> row(const long long& i) const;
    std::vector<T> col(const long long& j) const;

    // Resize the matrix.
    void resize(const long long& rs, const long long& cs, const T& def_val = T());

    // Accessors for dimensions and internal data.
    inline long long row_count() const;
    inline long long col_count() const;
    inline std::vector<T> data() const;
    inline std::vector<std::vector<T>> mat() const;
};

#include "../src/dmatreal.hpp"
