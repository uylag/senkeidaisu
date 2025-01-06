#include <iostream>
#include "../include/matrix.hpp"

int main() {
    std::vector<std::vector<long double>> matrix({
        {1.0L, 2.0L, 3.0L},
        {0.0L, 1.0L, 2.0L},
        {3.0L, 0.0L, 0.0L}
    });
    // SparseMatrix<long double> sparse;
    // SparseMatrix<long double>::CompressedColumnStorage crs = sparse.build_ccs(matrix);
    
    // std::cout << crs.get(1, 0);
    // crs(1, 0) = 5;
    // std::cout << crs.get(1, 0);
    Matrix<long double> mat(matrix);
}