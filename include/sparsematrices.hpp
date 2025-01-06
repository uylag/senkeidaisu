#pragma once
#include "matrix.hpp"
#include "densematrices.hpp"

template <typename T>
class SparseMatrix {
public:
    struct CompressedRowStorage {
        std::vector<T> values;
        std::vector<long long> col_ind;
        std::vector<long long> row_ptr;

        void set(const T& value, const long long& i, const long long& j);
        T& get(const long long& i, const long long& j);
        T& operator()(const long long& i, const long long& j);
        const T& operator()(const long long& i, const long long& j) const;
    };

    struct CompressedColumnStorage {
        std::vector<T> values;
        std::vector<long long> row_ind;
        std::vector<long long> col_ptr;

        void set(const T& value, const long long& i, const long long& j);
        T& get(const long long& i, const long long& j);
        T& operator()(const long long& i, const long long& j);
        const T& operator()(const long long& i, const long long& j) const;
    };
    
    SparseMatrix() = default;
    SparseMatrix(const std::vector<T>& c_matrix, const long long& i, const long long& j);
    SparseMatrix(const std::vector<std::vector<T>>& matrix);

    CompressedRowStorage build_crs(const std::vector<T>& c_matrix, const long long& i, const long long& j);
    CompressedRowStorage build_crs(const std::vector<std::vector<T>>& matrix);
    CompressedRowStorage build_crs(const DenseMatrix<T>& matrix);
    // CompressedRowStorage build_crs(Matrix<T>& matrix);

    CompressedColumnStorage build_ccs(const std::vector<T>& c_matrix, const long long& i, const long long& j);
    CompressedColumnStorage build_ccs(const std::vector<std::vector<T>>& matrix);
    CompressedColumnStorage build_ccs(const DenseMatrix<T>& matrix);
    // CompressedColumnStorage build_ccs(Matrix<T>& matrix);

protected:

private:
    CompressedRowStorage crs;
    CompressedColumnStorage ccs;
};

#include "../src/smatreal.hpp"