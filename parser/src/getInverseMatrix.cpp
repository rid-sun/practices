#include "getInverseMatrix.h"

/*
 * 采用LU分解法求得逆矩阵
 * @param JAC 所求逆矩阵的原矩阵
 * @param datum 规避节点编号
 */
void LU_decomposition(std::vector<std::vector<double>> &JAC)
{
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

    for (int i = 0; i < shape; i++)
    {
        W[i].resize(shape);
        W_n[i].resize(shape);
        L[i].resize(shape);
        L_n[i].resize(shape);
        U[i].resize(shape);
        U_n[i].resize(shape);
    }

    // 赋初值
    for (i = 0; i < shape; i++)
    {
        for (j = 0; j < shape; j++)
        {
            W[i][j] = JAC[i][j];
            L[i][j] = 0;
            U[i][j] = 0;
            L_n[i][j] = 0;
            U_n[i][j] = 0;
            W_n[i][j] = 0;
        }
    }

    // L对角置1
    for (i = 0; i < shape; i++)
    {
        L[i][i] = 1.0;
    }

    for (j = 0; j < shape; j++)
    {
        U[0][j] = W[0][j];
    }

    for (i = 1; i < shape; i++)
    {
        if (U[0][0] != 0)
            L[i][0] = W[i][0] / U[0][0];
    }

    for (i = 1; i < shape; i++)
    {
        // 求U
        for (j = i; j < shape; j++)
        {
            s = 0;
            for (k = 0; k < i; k++)
            {
                s += L[i][k] * U[k][j];
            }
            U[i][j] = W[i][j] - s;
        }
        // 求L
        for (d = i; d < shape; d++)
        {
            s = 0;
            for (k = 0; k < i; k++)
            {
                s += L[d][k] * U[k][i];
            }
            if (U[i][j] != 0)
                L[d][i] = (W[d][i] - s) / U[i][i];
        }
    }

    // 求L的逆
    for (j = 0; j < shape; j++)
    {
        for (i = j; i < shape; i++)
        {
            if (i == j)
            {
                if (L[i][j] != 0)
                    L_n[i][j] = 1 / L[i][j];
            }
            else if (i < j)
                L_n[i][j] = 0;
            else
            {
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
    for (i = 0; i < shape; i++)
    {
        for (j = i; j >= 0; j--)
        {
            if (i == j)
            {
                if (U[j][i] != 0)
                {
                    U_n[j][i] = 1 / U[j][i];
                }
            }
            else if (j > i)
                U_n[j][i] = 0;
            else
            {
                s = 0.;
                for (k = j + 1; k <= i; k++)
                {
                    s += U[j][k] * U_n[k][i];
                }
                if (U[j][j] != 0)
                {
                    U_n[j][i] = -1 / U[j][j] * s;
                }
            }
        }
    }

    for (i = 0; i < shape; i++)
    {
        for (j = 0; j < shape; j++)
        {
            for (k = 0; k < shape; k++)
            {
                W_n[i][j] += U_n[i][k] * L_n[k][j];
            }
        }
    }

    for (i = 0; i < shape; i++)
    {
        for (j = 0; j < shape; j++)
        {
            JAC[i][j] = W_n[i][j];
        }
    }
}

//得到给定矩阵src的逆矩阵保存到des中。
bool GetMatrixInverse(std::vector<std::vector<double>> &JAC)
{
    int n = JAC[0].size();
    double flag = getA(JAC, n);
    std::vector<std::vector<double>> t(n, std::vector<double>(n, 0));
    if (flag == 0)
    {
        return false;
    }
    else
    {
        getAStart(JAC, n, t);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                JAC[i][j] = t[i][j] / flag;
            }
        }
    }

    return true;
}

//按第一行展开计算|A|
double getA(std::vector<std::vector<double>> &JAC, int n) {
    if (n == 1) {
        return JAC[0][0];
    }
    double ans = 0;
    std::vector<std::vector<double>> temp(JAC[0].size(), std::vector<double>(JAC[0].size(), 0));
    int i, j, k;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n - 1; j++)
        {
            for (k = 0; k < n - 1; k++)
            {
                temp[j][k] = JAC[j + 1][(k >= i) ? k + 1 : k];
            }
        }
        double t = getA(temp, n - 1);
        if (i % 2 == 0)
        {
            ans += JAC[0][i] * t;
        }
        else
        {
            ans -= JAC[0][i] * t;
        }
    }
    return ans;
}
//计算每一行每一列的每个元素所对应的余子式，组成A*
void getAStart(std::vector<std::vector<double>> &JAC, int n, std::vector<std::vector<double>> &ans)
{
    if (n == 1)
    {
        ans[0][0] = 1;
        return;
    }
    int i, j, k, t;
    std::vector<std::vector<double>> temp(JAC[0].size(), std::vector<double>(JAC[0].size(), 0));
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            for (k = 0; k < n - 1; k++)
            {
                for (t = 0; t < n - 1; t++)
                {
                    temp[k][t] = JAC[k >= i ? k + 1 : k][t >= j ? t + 1 : t];
                }
            }

            ans[j][i] = getA(temp, n - 1);
            if ((i + j) % 2 == 1)
            {
                ans[j][i] = -ans[j][i];
            }
        }
    }
}