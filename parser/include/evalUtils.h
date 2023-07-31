#include "parser_v2.h"
#include "getInverseMatrix.h"
#include <Eigen/Dense>
#include <random>
#include <Python.h>

namespace evaluation{ 
    // 绘图相关
    vector<double> x_axis_data;
    vector<vector<double>> y_axis_data;
    
    NodeHead nodeList; // 注意，这里并不是引用类型，所以在初始化值的时候，他其实是网表对应内容的一块“复制品”
    CompHead compList;  // 换言之，如果改变这里的对象内容，并不会对网表实体中的内容造成任何的影响
    ModelHead modelList;
    // Netlist netlist;

    // 分析求解相关
    vector<double> lamda;
    vector<double> F_X, X, X_n;
    vector<double> G, a; // 同伦法相关使用的变量
    vector<vector<double>> JAC;
#ifdef USE_EIGEN_SJT
    Eigen::MatrixXd _JAC, _G;
    Eigen::VectorXd _F_X, _X, _X_n, _a;
#endif

    string outFileName, title;
    int datum, lastnode, step, total;
    double tran_stop, tran_initialVal;
    double tran_step;
    AnalysisType analysisType;

    const int ITERATIONNUMS = 15;
    const int RANDOM_MIN = 333333, RANDOM_MAX = 999999;
    const double ERRORGAP = 0.0001;
    Boolean isSuccess = FALSE;

    int plot(vector<double> X, vector<vector<double>> Y, string x_name, string y_name, string title, int step);
    void test_aq(vector<double>& X, vector<double>& a);
    void tranProcess();
    void newtonRaphson();
    void newtonIterHomo();
    void analysisProcess();
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
    tran_stop = netlist.getTranStop();
    tran_step = 1e-3;
    analysisType = netlist.getAnalysisType();
    title = netlist.getTitle();

    // 预处理，初始化方程矩阵、雅克比矩阵等
    X.resize(total - 1);
    F_X.resize(total - 1);
    X_n.resize(total - 1);
    JAC.resize(total - 1);
    G.resize(total - 1);
    a.resize(total - 1);
#ifdef USE_EIGEN_SJT
    _JAC.resize(total - 1, total - 1);
    _G.resize(total - 1, total - 1);
    _F_X.resize(total - 1);
    _X.resize(total - 1);
    _X_n.resize(total - 1);
    _a.resize(total - 1);
#endif
    

    random_device rd;
    default_random_engine engine(rd());
    uniform_int_distribution<int> distr(RANDOM_MIN, RANDOM_MAX);

    for (int i = 0; i < total - 1;i++) {
        X[i] = distr(engine) * 1.0 / 1000000; // 给X变量赋初值，随机选取初值
        a[i] = X[i];
    }
    fill(F_X.begin(), F_X.end(), 0);
    for (int i = 0; i < total - 1;i++){
        JAC[i].resize(total - 1);
        fill(JAC[i].begin(), JAC[i].end(), 0);
    }
    fill(G.begin(), G.end(), 1e-3);

#ifdef USE_EIGEN_SJT
    for (int i = 0; i < total - 1;i++) {
        _a(i) = _X(i) = distr(engine) * 1.0 / 1000000; // 给X变量赋初值，随机选取初值
    }
    _F_X.setZero();
    _JAC.setZero();
    _X_n.setZero();
    _G.diagonal().setConstant(1e-3);
#endif

    lamda.resize(11);
    for (int i = 0; i <= 10;i++) {
        lamda[i] = i * 1.0 / 10;
    }
}

// 牛顿迭代过程
void evaluation::newtonRaphson() {

    // 测试，赋予一定的初值
    test_aq(X, a);
    
    // 注意电压源引入补偿方程的初值是确定的
    // 2023/7/30 add：而且还需明确，只有一端接地，才能直接赋值的
    Component *tempComp = compList.getComp(0);
    while(tempComp != NULL) {
        if(tempComp->getType() == VSource) {
            int x1 = tempComp->getConVal(1) > datum ? tempComp->getConVal(1) - 1 : tempComp->getConVal(1);
            int x0 = tempComp->getConVal(0) > datum ? tempComp->getConVal(0) - 1 : tempComp->getConVal(0);
            if(tempComp->getConVal(0) == datum) {
                X[x1] = tempComp->getVal();
            }   
            else if(tempComp->getConVal(1) == datum) {
                X[x0] = tempComp->getVal();
            }
        }
        tempComp=tempComp->getNext();
    }

    // for (int i = 0; i < X.size(); i++) {
    //     cout << X[i] << " ";
    // }
    // cout << endl;

    // 2. 牛顿迭代开始
    for (int i = 0; i < ITERATIONNUMS; i++) {
        cout << endl
             << "     the beginning of " << i + 1 << "th iteration=========================" << endl;

        // 2.1 求出jacobian、F（X）等
        generateMatrix(nodeList, compList, modelList, F_X, X, JAC, outFileName, datum, lastnode, i + 1);

        // for (int i = 0; i < total - 1; i++) {
        //     for (int j = 0; j < total - 1; j++)
        //         cout << "J(" << i << "," << j << ") = " << JAC[i][j] << " ";
        //     cout << endl;
        // }
    
    #ifndef USE_EIGEN_SJT
        // 2.2 求出jacobian矩阵的逆矩阵
        LU_decomposition(JAC);

        cout << endl
             << "==========================================================================" << endl;

        // 2.3 开始牛顿法迭代计算
        for (int x = 0; x < total - 1; x++) {
            double temp = 0;
            for (int y = 0; y < total - 1; y++) {
                temp += JAC[x][y] * F_X[y];
            }
            X_n[x] = X[x] - temp;
        }
        cout << endl;

        // 2.4 进行迭代退出条件的判断
        Boolean isCorrect = TRUE;
        for (int x = 0; x < total - 1; x++) {
            if (abs(X[x] - X_n[x]) > ERRORGAP) {
                isCorrect = FALSE;
            }
            X[x] = X_n[x];
        }
        if (isCorrect == TRUE) {
            isSuccess = TRUE;
            break;
        }
        step++;

        // 2.5 重置为0，进行下一次的迭代
        for (int j = 0; j < total - 1; j++) {
            fill(JAC[j].begin(), JAC[j].end(), 0);
            F_X[j] = 0;
        }
    #else
        for (size_t hh = 0; hh < F_X.size(); ++hh) {
            _F_X(hh) = F_X[hh];
            _X(hh) = X[hh];
        }
        for (size_t hh = 0; hh < JAC.size(); ++hh) {
            for (size_t kk = 0; kk < JAC[0].size(); ++kk) {
                _JAC(hh, kk) = JAC[hh][kk];
            }
        }
        // _F_X = Eigen::Map<Eigen::VectorXd>(F_X.data(), F_X.size());
        // _X = Eigen::Map<Eigen::VectorXd>(X.data(), X.size());
        // _JAC = Eigen::MatrixXd::Map(JAC[0].data(), JAC.size(), JAC[0].size());
        _X_n = _X - _JAC.colPivHouseholderQr().solve(_F_X);

        if ((_X_n - _X).norm() < ERRORGAP) {
            isSuccess = TRUE;
            break;
        }

        step++;

        for (int j = 0; j < total - 1; j++) {
            fill(JAC[j].begin(), JAC[j].end(), 0);
            F_X[j] = 0;
            X[j] = _X_n(j);
        }
    #endif
    }

    if(isSuccess == TRUE) {
        cout << endl
             << "============In the " << step << "th iteration, the solution is successful==============" << endl;
        cout << endl
             << "========================vector x is followings: ==========================" << endl;
        for (int i = 0; i < total - 1; i++){
            cout << "X[" << (i >= datum ? i + 1 : i) << "] = " << X[i] << "      ";
        }
    } else {
        cout << endl
             << "==========Failed to resolve the problem within the specified maximum number of iterations===============" << endl;
    }
}

// 同伦法求解过程
void evaluation::newtonIterHomo() {

    // 测试，赋予一定的初值
    test_aq(X, a);

#ifdef USE_EIGEN_SJT
    _a = Eigen::Map<Eigen::VectorXd>(a.data(), a.size());
#endif

    // 注意电压源引入补偿方程的初值是确定的
    // 2023/7/30 add：而且还需明确，只有一端接地，才能直接赋值的
    Component *tempComp = compList.getComp(0);
    while(tempComp != NULL) {
        if(tempComp->getType() == VSource) {
            // if(tempComp->getConVal(0) != datum) {// 这里的设定有点糙
            //     X[tempComp->getConVal(0)] = tempComp->getVal();
            //     a[tempComp->getConVal(0)] = tempComp->getVal();
            // }   
            // else {
            //     X[tempComp->getConVal(1)] = tempComp->getVal();
            //     a[tempComp->getConVal(1)] = tempComp->getVal();
            // }
            int x1 = tempComp->getConVal(1) > datum ? tempComp->getConVal(1) - 1 : tempComp->getConVal(1);
            int x0 = tempComp->getConVal(0) > datum ? tempComp->getConVal(0) - 1 : tempComp->getConVal(0);
            if(tempComp->getConVal(0) == datum) {
                X[x1] = tempComp->getVal();
                a[x1] = tempComp->getVal();
            }   
            else if(tempComp->getConVal(1) == datum) {
                X[x0] = tempComp->getVal();
                a[x0] = tempComp->getVal();
            }
                
        }
        tempComp=tempComp->getNext();
    }
    
    // for (int i = 0; i < X.size(); i++) {
    //     cout << X[i] << " ";
    // }

    // 对要绘图的量进行清空
    x_axis_data.clear();
    y_axis_data.clear();
    y_axis_data.resize(total - 1);

    // 开始同伦法求解迭代
    for (int i = 1; i <= 10;i++) {
        cout << endl
             << "lamda = " << lamda[i] << " iteration---" << endl;
        step = 1;
        while(step < ITERATIONNUMS) {
            cout << endl
                 << "     step " << step << "th iteration ---" << endl;
            Boolean isSuccess = TRUE;
            
            // 1. 计算jacobian矩阵 和 F_X
            generateMatrix(nodeList, compList, modelList, F_X, X, JAC, outFileName, datum, lastnode, step);
            // for (int i = 0; i < total - 1; i++) {
            //     for (int j = 0; j < total - 1; j++)
            //         cout << "J(" << i << "," << j << ") = " << JAC[i][j] << " ";
            //     cout << endl;
            // }
        
        #ifndef USE_EIGEN_SJT
            // 2. 计算H(x, λ) 【这里可能需要再前面得先进行一次求解FX】
            // H(x, λ) = (1 - λ)G(x - a) + λF(x)
            for (int x = 0; x < total - 1; x++) {
                X_n[x] = (1 - lamda[i]) * (1e-3) * (X[x] - a[x]) + lamda[i] * F_X[x];
                for (int y = 0; y < total - 1; y++) {
                    JAC[x][y] = lamda[i] * JAC[x][y] + (x == y ? (1 - lamda[i]) * 1e-3 : 0);
                }
                if(abs(X_n[x]) > ERRORGAP) { // 注意这个误差的处理
                    isSuccess = FALSE;
                }
            }
            if(isSuccess == TRUE) break;
            
            // 3. 迭代下一个解
            // 3.1 JAC矩阵求逆
            LU_decomposition(JAC);

            // 3.2 求解下一个解
            vector<double> temp(X_n);
            for (int x = 0; x < total - 1; x++) {
                double t = 0;
                for (int y = 0; y < total - 1; y++) {
                    t += JAC[x][y] * temp[y];
                }
                X_n[x] = X[x] - t;
            }

            // 4. 重置 and 赋值
            for (int h = 0; h < total - 1; h++) {
                X[h] = X_n[h];
                fill(JAC[h].begin(), JAC[h].end(), 0);
                F_X[h] = 0;
            }
        #else
            for (size_t hh = 0; hh < F_X.size(); ++hh) {
                _F_X(hh) = F_X[hh];
                _X(hh) = X[hh];
            }
            // _F_X = Eigen::Map<Eigen::VectorXd>(F_X.data(), F_X.size());
            // _X = Eigen::Map<Eigen::VectorXd>(X.data(), X.size());
            // _JAC = Eigen::Map<Eigen::MatrixXd>(&(JAC[0][0]), JAC.size(), JAC[0].size()); // 这种赋值方式感觉会出错，因为地址是不连续的码？
            for (size_t hh = 0; hh < JAC.size(); ++hh) {
                for (size_t kk = 0; kk < JAC[0].size(); ++kk) {
                    _JAC(hh, kk) = lamda[i] * JAC[hh][kk] + (hh == kk ? (1 - lamda[i]) * 1e-3 : 0);
                }
            }
            // cout << _F_X << "----" << endl << _X << "----" << endl << _JAC << endl;

            // _X = _X - _JAC.colPivHouseholderQr().solve(_F_X);

            _X_n = (1 - lamda[i]) * _G * (_X - _a) + lamda[i] * _F_X;

            if (_X_n.norm() < ERRORGAP) break;

            _JAC = _JAC.inverse();

            _X = _X - _JAC * _X_n;

            for (int h = 0; h < total - 1; h++) {
                X[h] = _X(h);
                fill(JAC[h].begin(), JAC[h].end(), 0);
                F_X[h] = 0;
            }
        #endif
            cout << endl
                 << "     step " << step << "th iteration end.----" << endl;
            step++;
        }
        cout << endl
             << "================the end of lamda = " << lamda[i] << " iteration=============================" << endl;
        
        // 设置要绘图的量
        x_axis_data.push_back(lamda[i]);
        for (int j = 0, h = 0; j < total - 1; j++) {
            y_axis_data[h++].push_back(X[j]);
        }
    }

    cout << endl
         << "========================vector x is followings: ==========================" << endl;
    for (int i = 0; i < total - 1; i++){
        cout << "X[" << (i >= datum ? i + 1 : i) << "] = " << X[i] << "      ";
    }

    // 开始绘图
    // 绘图
    if (plot(x_axis_data, y_axis_data, "lamda", "voltage(v)", title, 0) != 512) {
        cerr << endl
             << "Can't draw the fig!" << endl;
    }
}

// 进行tran分析的处理
/* TODO：需要做一个更普适性的工作 */
/* 目前有以下限制：
/* 1. 激励信号默认是从0开始，即初值为0             */
/* 2. 步长没有动态调整，固定为1ms                 */
/* 3. 在电路中只包含电容、直流电源、电阻，且电容和电源数目为1个，并且要求一端接地 */
/* 4. 对于误差的计算也没有做处理                  */
void evaluation::tranProcess() {
    cout << endl
         << "=========================Tran process has started!========================"<<endl;
    double stop = tran_stop, ans = 0, C_value = 0;
    double step = tran_step, total_resistor1 = 0, total_resistor2 = 0;
    int id1 = NA, id2 = NA;
    Component *comp = compList.getComp(0);
    while (comp != NULL) {
        if (comp->getType() == Capacitor) {
            id1 = comp->getConVal(0);
            id2 = comp->getConVal(1);
            break;
        }
        comp = comp->getNext();
    }
    if (id1 == NA || id2 == NA) {
        cerr << "There is no capacitor, so you can't do transient analysis!";
        exit(1);
    }
    id1 = id1 == datum ? id2 : id1;
    Node *node = nodeList.getNode(0);
    while (node != NULL) {
        if (node->getNameNum() != datum) {
            Connections *conList =  node->getConList();
            while (conList != NULL) {
                if (conList->comp->getType() == Resistor && conList->comp->getConVal(0) == node->getNameNum()) {
                    if (node->getNameNum() == id1)
                        total_resistor1 += 1 / conList->comp->getVal();
                    else
                        total_resistor2 += 1 / conList->comp->getVal();
                } else if (conList->comp->getType() == VSource) {
                    tran_initialVal = conList->comp->getVal();
                } else if (conList->comp->getType() == Capacitor) {
                    C_value = conList->comp->getVal();
                }
                conList = conList->next;
            }
        }
        node = node->getNext();
    }
    x_axis_data.clear();
    y_axis_data.clear();
    y_axis_data.push_back({});
    for (int h = 1; h <= (int)(stop / step); h++) {
        // ans = (tran_initialVal * total_resistor2 + C_value * ans / step) / (C_value / step + 1 / total_resistor1 + 1);
        ans = (tran_initialVal + C_value * ans / (step * total_resistor2)) / (C_value / (step * total_resistor2) + 1 + total_resistor1 / total_resistor2);
        x_axis_data.push_back(h * step);
        y_axis_data.back().push_back(ans);
    }
    cout << endl
         << "The value at that time is " << ans << endl;
    
    // 绘图
    if (plot(x_axis_data, y_axis_data, "time(s)", "voltage(v)", title, 0) != 512) {
        cerr << endl
             << "Can't draw the fig!" << endl;
    }
}

// 分析过程的入口
void evaluation::analysisProcess() {
    switch (analysisType) {
    case TRAN:
        tranProcess();
        break;
    case DC:
        newtonRaphson();
        // newtonIterHomo();
        break;
    case AC:
        break;
    }
}

void evaluation::test_aq(vector<double> &X, vector<double> &a) {
    vector<vector<double>> all_inival;
    string path = "..\\testcase\\testcase3\\Initial_val3.txt"; // 文件路径
    char *line = (char *)malloc(sizeof(char) * 1024);
    char *p;
    FILE *fp = fopen(path.c_str(), "r");
    if (fp == NULL) {
        cout << "File opening error!" << endl;
        return;
    }
    while (!feof(fp)) {
        fgets(line, 1024, fp);
        p = strtok(line, " ");
        if (isdigit(p[0])) {
            vector<double> t;
            // t.push_back(0);
            while (p != NULL) {
                t.push_back(strtod(p, NULL));
                p  = strtok(NULL, " ");
            }
            all_inival.push_back(t);
        }
    }
    fclose(fp);
    
    // 选取所有初值中的一种来赋值
    for (int i = 0; i < all_inival[0].size(); i++) {
        a[i] = X[i] = all_inival[0][i];
        // cout << X[i] << " ";
    }
    // cout << X.size() << endl;

    return;
}

// 绘图
int evaluation::plot(vector<double> X, vector<vector<double>> Y, string x_name, string y_name, string title, int step) {
    
    // 很重要
    Py_SetPythonHome(L"D:/Anaconda3");
    
    // 初始化python环境
    Py_Initialize();
    if (!Py_IsInitialized()) {
        printf("Py_Initialize failed!!!\n");
        return -1;
    }

    // 测试语句
    // PyRun_SimpleString("print('Hello Python!')\n");
    PyRun_SimpleString("import os,sys");//执行import语句，把当前路径的上一级加入路径中，为了找到plot.py
    PyRun_SimpleString("sys.path.append('../scripts')");
    // PyRun_SimpleString("print(os.getcwd())");//测试打印当前路径

    // 定义调用函数时的相关变量
    PyObject *pModule;
    PyObject *pFunction;
    PyObject *pArgs;
    PyObject *pRetValue;

    // 加载脚本进来
    pModule = PyImport_ImportModule("plot"); // 注意这里一定要是plot，不带.py
    if (!pModule) {
        printf("import python failed!!!\n");
        return -1;
    }

    // 查找工具函数
    pFunction = PyObject_GetAttrString(pModule, "plot_util");
    if (!pFunction) {
        printf("get python function failed!!!\n");
        return -1;
    }

    // 构建参数
    // 预处理vector(cpp) -> list(py)
    PyObject *list1 = PyList_New(X.size()); // X -> x_axis_data
    for (int i = 0; i < X.size(); i++) {
        PyList_SetItem(list1, i, Py_BuildValue("d", X[i]));
    }
    PyObject *list2 = PyList_New(0);
    for (int i = 0; i < Y.size(); i++) {
        PyObject *t = PyList_New(Y[i].size());
        for (int j = 0; j < Y[i].size(); j++) {
            PyList_SetItem(t, j, Py_BuildValue("d", Y[i][j]));
        }
        PyList_Append(list2, t);
        Py_DECREF(t);   
    }
    
    pArgs = PyTuple_New(6);
    PyTuple_SetItem(pArgs, 0, list1);
    PyTuple_SetItem(pArgs, 1, list2);
    PyTuple_SetItem(pArgs, 2, Py_BuildValue("s", x_name.c_str()));
    PyTuple_SetItem(pArgs, 3, Py_BuildValue("s", y_name.c_str()));
    PyTuple_SetItem(pArgs, 4, Py_BuildValue("s", title.c_str()));
    PyTuple_SetItem(pArgs, 5, Py_BuildValue("s", ("fig" + to_string(step)).c_str()));

    // 调用函数
    pRetValue = PyObject_CallObject(pFunction, pArgs);

    // 清空PyObject
    Py_DECREF(pModule);
    Py_DECREF(pFunction);
    Py_DECREF(pArgs);
    Py_DECREF(pRetValue);
    Py_DECREF(list1);
    Py_DECREF(list2);

    // 终止Py环境
    Py_Finalize();

    return 512;
}