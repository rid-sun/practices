#include <iostream>
#include <vector>

void LU_decomposition(std::vector<std::vector<double>> &JAC);

bool GetMatrixInverse(std::vector<std::vector<double>> &JAC);

double getA(std::vector<std::vector<double>> &JAC, int n);
void getAStart(std::vector<std::vector<double>> &JAC, int n, std::vector<std::vector<double>> &ans);
