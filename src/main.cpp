#include <bits/stdc++.h>
#include "../include/matrix.hpp"
using namespace std;

void swap(int* a, int* b)
{
    int* temp = a;
    a = b;
    b = temp;
}

int main() {
    Matrix<long double> matrix({
        {1.0L, 2.0L, 3.0L},
        {0.0L, 1.0L, 2.0L},
        {3.0L, 0.0L, 0.0L}
    });
    
    matrix = matrix + matrix;

    std::cout << matrix(2, 0);
}