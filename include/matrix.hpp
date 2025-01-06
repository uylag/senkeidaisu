#pragma once
#include "densematrices.hpp"
#include "sparsematrices.hpp"

template <typename T>
class Matrix {
    public:
        Matrix(const std::vector<T>& c_matrix, const long long& i, const long long& j);
        Matrix(const std::vector<std::vector<T>> matrix);
    protected:

    private:
        SparseMatrix<T> sparse;
        DenseMatrix<T> dense;
};

#include "../src/matreal.hpp"