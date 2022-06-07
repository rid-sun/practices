#include "getInverseMatrix.h"

/*
 * 采用LU分解法求得逆矩阵
 * @param JAC 所求逆矩阵的原矩阵
 */
void LU_decomposition(std::vector<std::vector<double>> &JAC) {
    int shape = JAC[0].size();
    int i, j, k, d;
    double s;
    std::vector<std::vector<double>> W, W_n, L, U, L_n, U_n;
    W.resize(shape);
    W_n.resize(shape);
    L.resize(shape);
    L_n.resize(shape);
    U.resize(shape);
    U_n.resize(shape);

    for (int i = 0; i < shape;i++){
        W[i].resize(shape);
        W_n[i].resize(shape);
        L[i].resize(shape);
        L_n[i].resize(shape);
        U[i].resize(shape);
        U_n[i].resize(shape);
    }

    // 赋初值
    for (i = 0; i < shape; i++) {
        for (j = 0; j < shape; j++) {
            W[i][j] = JAC[i][j];
            L[i][j] = 0;
            U[i][j] = 0;
            L_n[i][j] = 0;
            U_n[i][j] = 0;
            W_n[i][j] = 0;
        }
    }

    // L对角置1
    for (i = 0; i < shape; i++) {
        L[i][i] = 1.0;
    }

    for (j = 0; j < shape; j++) {
        U[0][j] = W[0][j];
    }

    for (i = 1; i < shape; i++) {
        L[i][0] = W[i][0] / U[0][0];
    }

    for (i = 1; i < shape; i++) {
        // 求U
        for (j = i; j < shape; j++) {
            s = 0;
            for (k = 0; k < i; k++) {
                s += L[i][k] * U[k][j];
            }
            U[i][j] = W[i][j] - s;
        }
        // 求L
        for (d = i; d < shape; d++) {
            s = 0;
            for (k = 0; k < i; k++) {
                s += L[d][k] * U[k][i];
            }
            L[d][i] = (W[d][i] - s) / U[i][i];
        }
    }

    // 求L的逆
    for (j = 0; j < shape; j++) {
        for (i = j; i < shape; i++) {
            if (i == j)
                L_n[i][j] = 1 / L[i][j];
            else if (i < j)
                L_n[i][j] = 0;
            else {
                s = 0.;
                for (k = j; k < i; k++)
                {
                    s += L[i][k] * L_n[k][j];
                }
                L_n[i][j] = -L_n[j][j] * s;
            }
        }
    }

    //求U的逆
    for (i = 0; i < shape; i++) {
        for (j = i; j >= 0; j--) {
            if (i == j)
                U_n[j][i] = 1 / U[j][i];
            else if (j > i)
                U_n[j][i] = 0;
            else {
                s = 0.;
                for (k = j + 1; k <= i; k++)
                {
                    s += U[j][k] * U_n[k][i];
                }
                U_n[j][i] = -1 / U[j][j] * s;
            }
        }
    }

    for (i = 0; i < shape; i++) {
        for (j = 0; j < shape; j++) {
            for (k = 0; k < shape; k++) {
                W_n[i][j] += U_n[i][k] * L_n[k][j];
            }
        }
    }

    for (i = 0; i < shape; i++) {
        for (j = 0; j < shape; j++) {
            JAC[i][j] = W_n[i][j];
        }
    }
}