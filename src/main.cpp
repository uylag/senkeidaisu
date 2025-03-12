#include "bits/stdc++.h"
#include <senkeidaisu.hpp>
#include <vector>
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
    senkeidaisu::Matrix<long double> m1({
        {1, 2, 3},
        {4, 5, 6},
        {7, 8, 0}
    });

    std::cout << m1.det();

    return 0;
}