#include <iostream>
#include <sstream>
#include <vector>
#include <random>
#include "parser_v2.h"
#include "getInverseMatrix.h"
#include "config.h"
using namespace std;

// 分析求解相关
vector<double> F_X, X, X_n;
vector<vector<double>> JAC;
NodeHead nodeList;
CompHead compList;
ModelHead modelList;
int datum, lastnode, step, total;
string inFileName, outFileName;

const int ITERATIONNUMS = 15;
const int RANDOM_MIN = 333333, RANDOM_MAX = 999999;
const double ERRORGAP = 0.0001;
Boolean isSuccess = FALSE;

int main(int argc, char** argv) {
    config::parse_args(argc, argv);
    if(config::instrError) {
        cout << argv[0] << " Invalid instruction, please use '-h' to understand the usage!" << endl;
    } else if(config::help) {
        cout << "Usage: " << argv[0] << " -f FILENAME <-d datum> <-h> <-o OUTPUT_FILE>" << endl;
        cout << "  ------- Listing options -------" << endl;
        cout << "  -h              Print usage and this help message." << endl;
        cout << "  -o[- filename]  When is '-', it means the output stream is standard output; When 'filename', represents the output stream as the file." << std::endl;
        cout << "  -f[- filename]  When is '-', it means the input stream is standard output; When 'filename', represents the input stream as the file." << std::endl;
        cout << "  -d              Specify the start node number." << endl;
    } else {
        // cout<<"datum: "<<config::datum<<" inFileName: "<<config::inFileName<<" outFileName: "<<config::outFileName<<endl;
        // 参数、变量等初始化
        datum = config::datum;
        lastnode = NA;
        step = 1;
        inFileName = config::inFileName;
        outFileName = config::outFileName;

        // 解析网表文件，存储信息到nodeList、compList、modelList中
        parseNetList(nodeList, compList, modelList, datum, lastnode, inFileName, outFileName);

        /* 实则lastnode + 1 == nodeList.getCount() */
        total = lastnode + compList.getCount(VSource) + 1;
        // 预处理，初始化方程矩阵、雅克比矩阵等
        X.resize(total);
        F_X.resize(total);
        X_n.resize(total);
        JAC.resize(total);
        random_device rd;
        default_random_engine engine(rd());
        uniform_int_distribution<int> distr(RANDOM_MIN, RANDOM_MAX);
        for (int i = 0; i < total;i++){
            X[i] = distr(engine) * 1.0 / 1000000; // 给X变量赋初值，随机选取初值
        }
        fill(F_X.begin(), F_X.end(), 0);
        for (int i = 0; i < total;i++){
            JAC[i].resize(total);
            fill(JAC[i].begin(), JAC[i].end(), 0);
        }
        
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
        // for(int i=0;i<X.size();i++){
        //     cout<<"X["<<i<<"] = "<<X[i]<<" ";
        // }
        // cout<<endl;
        // for(int i=0;i<X.size();i++){
        //     cout<<"F_X["<<i<<"] = "<<F_X[i]<<" ";
        // }
        // cout<<endl;
        // for(int i=0;i<X.size();i++){
        //     for(int j=0;j<X.size();j++){
        //         cout<<"JAC["<<i<<","<<j<<"] = "<<JAC[i][j]<<" ";
        //     }
        //     cout<<endl;
        // }

        // 迭代开始
        for (int i = 0; i < ITERATIONNUMS; i++) {
            cout << "===================the beginning of " << step << "th iteration=========================" << endl;
            
            // for (int i = 0; i < total;i++){
            //     if(i == datum) continue;
            //     cout << "X[" << i << "] = " << X[i] << "      ";
            // }
            // cout << endl;

            generateMatrix(nodeList, compList, modelList, F_X, X, JAC, outFileName, datum, lastnode, step);

            // for (int i = 0; i < total; i++) {
            //     for (int j = 0; j < total; j++)
            //         cout << "J(" << i << "," << j << ") = " << JAC[i][j] << " ";
            //     cout << endl;
            // }

            vector<vector<double>> tJAC(total - 1, vector<double>(total - 1, 0));
            for (int u = 0, x = 0; u < total; u++) {
                if(u == datum) continue;
                for (int w = 0, y = 0; w < total;w++) {
                    if(w == datum) continue;
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

            // for (int i = 0; i < total;i++){
            //     if(i == datum) continue;
            //     cout << "X[" << i << "] = " << X[i] << "      ";
            // }
            // cout << endl;

            // 开始牛顿法迭代计算
            for (int h = 0, x = 0; h < total;h++) {
                if(h == datum) continue;
                double temp = 0;
                for (int k = 0, y = 0; k < total;k++){
                    if(k == datum) continue;
                    temp += tJAC[x][y] * F_X[k];
                    y++;
                }
                X_n[h] = X[h] - temp;
                x++;
            }
            cout << endl;

            Boolean isCorrect = TRUE;
            for (int j = 0; j < total;j++) {
                if(j == datum) continue;
                if(abs(X[j] - X_n[j]) > ERRORGAP){
                    isCorrect = FALSE;
                }
                X[j] = X_n[j];
            }
            if(isCorrect == TRUE) {
                isSuccess = TRUE;
                break;
            }
            step++;
            // 重置为0
            for (int j = 0; j < total;j++){
                fill(JAC[j].begin(), JAC[j].end(), 0);
                F_X[j] = 0;
            }
        }
        if(isSuccess == TRUE) {
            cout << "=======In the "<< step <<"th iteration, the solution is successful==============" << endl;
            cout << "=============vector x is followings: ==============" << endl;
            for (int i = 0; i < total;i++){
                if(i == datum) continue;
                cout << "X[" << i << "] = " << X[i] << "      ";
            }
        } else {
            cout << "==========Failed to resolve the problem within the specified maximum number of iterations===============" << endl;
        }
    }
    return 0;
}