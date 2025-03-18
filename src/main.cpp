// #include "bits/stdc++.h"
#include <senkeidaisu.hpp>
#include <lin_reg.hpp>
#include <vector>
#include <iostream>
using namespace std;

template <typename TN>
void cout_vector(std::vector<TN> vector)
{
    for (auto& elem: vector)
    {
        std::cout << elem << " ";
    }
    std::cout << '\n';
}

int main() 
{
    senkeidaisu::Matrix<long double> X({
        {1, 3, 5},
        {2, 4, 8},
        {3, 6, 9}
    });

    senkeidaisu::Matrix<long double> y({
        {6, 10, 12}
    });

    LinearRegression<long double> lr(X, y, 0.01);
    lr.fit(X, y);
    cout << lr.predict({4, 8, 10});

    return 0;
}