#include "parser_v2.h"
#include "getInverseMatrix.h"
#include <random>

namespace evaluation{
    // 分析求解相关
    vector<double> F_X, X, X_n, G, lamda, a;
    vector<vector<double>> JAC;
    NodeHead nodeList;
    CompHead compList;
    ModelHead modelList;
    int datum, lastnode, step, total;
    string outFileName;

    const int ITERATIONNUMS = 15;
    const int RANDOM_MIN = 333333, RANDOM_MAX = 999999;
    const double ERRORGAP = 0.0001;
    Boolean isSuccess = FALSE;

    void newtonRaphson();
    void newtonIterHomo();
    void preparation(Netlist &netlist, string outFileName);
}

// 初始化命名空间内的变量及做一些预处理
void evaluation::preparation(Netlist &netlist, string outFileName) {
    nodeList = netlist.getNodeHead();
    compList = netlist.getCompHead();
    modelList = netlist.getModelHead();
    datum = netlist.getDatum();
    lastnode = netlist.getLastnode();
    total = lastnode + compList.getCount(VSource) + 1; /* 实则lastnode + 1 == nodeList.getCount() */
    evaluation::outFileName = outFileName;

    // 预处理，初始化方程矩阵、雅克比矩阵等
    X.resize(total);
    F_X.resize(total);
    X_n.resize(total);
    JAC.resize(total);
    G.resize(total - 1);
    a.resize(total);
    random_device rd;
    default_random_engine engine(rd());
    uniform_int_distribution<int> distr(RANDOM_MIN, RANDOM_MAX);
    for (int i = 0; i < total;i++){
        X[i] = distr(engine) * 1.0 / 1000000; // 给X变量赋初值，随机选取初值
        a[i] = X[i];
    }
    fill(F_X.begin(), F_X.end(), 0);
    for (int i = 0; i < total;i++){
        JAC[i].resize(total);
        fill(JAC[i].begin(), JAC[i].end(), 0);
    }
    fill(G.begin(), G.end(), 1e-3);
    lamda.resize(11);
    for (int i = 0; i <= 10;i++) {
        lamda[i] = i * 1.0 / 10;
    }
}

// 牛顿迭代过程
void evaluation::newtonRaphson() {
    
    // 1. 设置X初始值
    // 注意电压源引入补偿方程的初值是确定的
    Component *tempComp = compList.getComp(0);
    while(tempComp != NULL) {
        if(tempComp->getType() == VSource) {
            if(tempComp->getConVal(0) != datum) // 这里的设定有点糙
                X[tempComp->getConVal(0)] = tempComp->getVal();
            else
                X[tempComp->getConVal(1)] = tempComp->getVal();
        }
        tempComp=tempComp->getNext();
    }

    // // 打印测试
    // for (int i = 0; i < X.size(); i++) {
    //     cout << "X[" << i << "] = " << X[i] << " ";
    // }
    // cout << endl;
    // for (int i = 0; i < X.size(); i++) {
    //     cout << "F_X[" << i << "] = " << F_X[i] << " ";
    // }
    // cout << endl;
    // for (int i = 0; i < X.size(); i++) {
    //     for (int j = 0; j < X.size(); j++) {
    //         cout << "JAC[" << i << "," << j << "] = " << JAC[i][j] << " ";
    //     }
    //     cout << endl;
    // }

    // 2. 牛顿迭代开始
    for (int i = 0; i < ITERATIONNUMS; i++) {
        cout << endl
             << "===================the beginning of " << i + 1 << "th iteration=========================" << endl;

        // for (int i = 0; i < total;i++) {
        //     if(i == datum) continue;
        //     cout << "X[" << i << "] = " << X[i] << "      ";
        // }
        // cout << endl;

        // 2.1 求出jacobian、F（X）等
        generateMatrix(nodeList, compList, modelList, F_X, X, JAC, outFileName, datum, lastnode, i + 1);

        // for (int i = 0; i < total; i++) {
        //     for (int j = 0; j < total; j++)
        //         cout << "J(" << i << "," << j << ") = " << JAC[i][j] << " ";
        //     cout << endl;
        // }

        // 2.2 求出jacobian矩阵的逆矩阵
        vector<vector<double>> tJAC(total - 1, vector<double>(total - 1, 0));
        for (int u = 0, x = 0; u < total; u++) {
            if (u == datum) continue;
            for (int w = 0, y = 0; w < total; w++) {
                if (w == datum) continue;
                tJAC[x][y] = JAC[u][w];
                y++;
            }
            x++;
        }

        // for (int i = 0; i < total - 1; i++) {
        //     for (int j = 0; j < total - 1; j++)
        //         cout << "tJ(" << i << "," << j << ") = " << tJAC[i][j] << " ";
        //     cout << endl;
        // }

        // GetMatrixInverse(tJAC);
        LU_decomposition(tJAC);

        cout << endl
             << "==========================================================================" << endl;

        // for (int i = 0; i < total - 1; i++) {
        //     for (int j = 0; j < total - 1; j++)
        //         cout << "tJ(" << i << "," << j << ") = " << tJAC[i][j] << " ";
        //     cout << endl;
        // }

        // for (int i = 0; i < total; i++) {
        //     if(i == datum) continue;
        //     cout << "X[" << i << "] = " << X[i] << "      ";
        // }
        // cout << endl;

        // 2.3 开始牛顿法迭代计算
        for (int h = 0, x = 0; h < total; h++) {
            if (h == datum) continue;
            double temp = 0;
            for (int k = 0, y = 0; k < total; k++) {
                if (k == datum) continue;
                temp += tJAC[x][y] * F_X[k];
                y++;
            }
            X_n[h] = X[h] - temp;
            x++;
        }
        cout << endl;

        // 2.4 进行迭代退出条件的判断
        Boolean isCorrect = TRUE;
        for (int j = 0; j < total; j++) {
            if (j == datum) continue;
            if (abs(X[j] - X_n[j]) > ERRORGAP) {
                isCorrect = FALSE;
            }
            X[j] = X_n[j];
        }
        if (isCorrect == TRUE) {
            isSuccess = TRUE;
            break;
        }
        step++;

        // 2.5 重置为0，进行下一次的迭代
        for (int j = 0; j < total; j++) {
            fill(JAC[j].begin(), JAC[j].end(), 0);
            F_X[j] = 0;
        }
    }

    if(isSuccess == TRUE) {
        cout << endl
             << "============In the " << step << "th iteration, the solution is successful==============" << endl;
        cout << endl
             << "========================vector x is followings: ==========================" << endl;
        for (int i = 0; i < total;i++){
            if(i == datum) continue;
            cout << "X[" << i << "] = " << X[i] << "      ";
        }
    } else {
        cout << endl
             << "==========Failed to resolve the problem within the specified maximum number of iterations===============" << endl;
    }
}

// 同伦法求解过程
void evaluation::newtonIterHomo() {

    // 注意电压源引入补偿方程的初值是确定的
    Component *tempComp = compList.getComp(0);
    while(tempComp != NULL) {
        if(tempComp->getType() == VSource) {
            if(tempComp->getConVal(0) != datum) {// 这里的设定有点糙
                X[tempComp->getConVal(0)] = tempComp->getVal();
                a[tempComp->getConVal(0)] = tempComp->getVal();
            }   
            else {
                X[tempComp->getConVal(1)] = tempComp->getVal();
                a[tempComp->getConVal(1)] = tempComp->getVal();
            }
                
        }
        tempComp=tempComp->getNext();
    }

    // 开始同伦法求解迭代
    for (int i = 1; i <= 10;i++) {
        cout << endl
             << "================the beginning of lamda = " << lamda[i] << "iteration======================" << endl;
        step = 1;
        while(step < ITERATIONNUMS) {
            cout << endl
                 << "============the beginning of " << step << "th iteration when lamda = " << lamda[i] << "================" << endl;
            Boolean isSuccess = TRUE;
            
            // 1. 计算jacobian矩阵 和 F_X
            generateMatrix(nodeList, compList, modelList, F_X, X, JAC, outFileName, datum, lastnode, step);
            
            // 2. 计算H(x, λ) 【这里可能需要再前面得先进行一次求解FX】
            // H(x, λ) = (1 - λ)G(x - a) + λF(x)
            for (int j = 0; j < total;j++) {
                if(j == datum) continue;
                X_n[j] = (1 - lamda[i]) * (1e-3) * (X[j] - a[j]) + lamda[i] * F_X[j];
                if(abs(X_n[j]) > ERRORGAP) { // 注意这个误差的处理
                    isSuccess = FALSE;
                }
            }
            if(isSuccess == TRUE) break;
            
            // 3. 迭代下一个解
            // 3.1 tJAC矩阵求逆
            vector<vector<double>> tJAC(total - 1, vector<double>(total - 1, 0));
            for (int u = 0, x = 0; u < total; u++) {
                if (u == datum) continue;
                for (int w = 0, y = 0; w < total; w++) {
                    if (w == datum) continue;
                    tJAC[x][y] = lamda[i] * JAC[u][w] + (u == w ? (1 - lamda[i]) * 1e-3 : 0);
                    y++;
                }
                x++;
            }
            LU_decomposition(tJAC);

            // 3.2 求解下一个解
            vector<double> temp(X_n);
            for (int h = 0, x = 0; h < total; h++) {
                if (h == datum) continue;
                double t = 0;
                for (int k = 0, y = 0; k < total; k++) {
                    if (k == datum) continue;
                    t += tJAC[x][y] * temp[k];
                    y++;
                }
                X_n[h] = X[h] - t;
                x++;
            }

            // 4. 重置 and 赋值
            for (int h = 0; h < total; h++) {
                X[h] = X_n[h];
                fill(JAC[h].begin(), JAC[h].end(), 0);
                F_X[h] = 0;
            }
            cout << endl
                 << "=================the end of " << step << "th iteration when lamda = " << lamda[i] << "==================" << endl;
            step++;
        }
        cout << endl
             << "================the end of lamda = " << lamda[i] << "iteration=============================" << endl;
    }

    cout << endl
         << "========================vector x is followings: ==========================" << endl;
    for (int i = 0; i < total; i++) {
        if (i == datum) continue;
        cout << "X[" << i << "] = " << X[i] << "      ";
    }
}