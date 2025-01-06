#pragma once
#include <vector>

template <typename T>
class DenseMatrix {
public:
    DenseMatrix() = default;
    DenseMatrix(const std::vector<T>& c_matrix, const long long& i, const long long& j);
    DenseMatrix(const std::vector<std::vector<T>>& matrix);

    std::vector<T> row(long long i) const;
    std::vector<T> col(long long j) const;

    inline long long row_count() const;
    inline long long col_count() const;
    inline std::vector<T> mat() const;

protected:

private:
    long long rows, columns;
    std::vector<T> matrix;
};

#include "../src/dmatreal.hpp"