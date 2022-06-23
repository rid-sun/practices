//  PARSER.H
//  Header for Parser program.

#ifndef PARSER_V2
#define PARSER_V2

#include <ctype.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <unordered_map>
using namespace std;

// 常量声明
const double K = 1.38E-23;
const double Q = 1.60E-19;
const double N = 38.78;
const int NameLength = 80, BufLength = 3000, NA = -1;

enum CompType {// 器件类型
    MOSFET,
    BJT,
    VSource,
    ISource,
    Inductor,
    Resistor,
    Diode,
    Capacitor
}; 
enum TranType {// tran类型
    NMOS,
    PMOS,
    NPN,
    PNP
};
enum Flag {// 标志【用于设置器件端口状态】
    UNSET,
    SET
};
enum Boolean {// 逻辑判断
    FALSE,
    TRUE
};
enum EquaType {// 方程类型
    Nodal,
    Modified
};
enum AnalysisType {// 分析类型
    DC,
    AC,
    TRAN
};

// 七种封装类
class Component;
class ComponentHead;
class Node;
class NodeHead;
class Model;
class ModelHead;
class Netlist;

// 两种信息类
struct Connectors {// 器件“端口”类封装信息
    Flag flag; // “端口”状态
    Node *node; // 本器件该“端口”所连节点实体
    int conNum; // 本器件该“端口”所连节点编号
};
struct Connections {// 连接关系信息封装体
    Connections *next; // 下一个“连接关系”
    Component *comp; // 指向器件实体
    int conNum; // 该关系所连节点的“端口”编号
};

class Component {// 器件
public:
    Component(CompType typeIn, double valueIn, double tempIn, int con0In,
              int con1In, int con2In, int con3In, Model *modelIn, char *nameIn);
    ~Component();
    CompType getType();
    Component *getNext();
    int getcompNum();
    int getNum();
    int getNodeNum(int conNum);
    double getVal();
    Node *getNode(int conNum);
    char *getName();
    int getConVal(int conNum);

    void setNext(Component *nextIn);void setNum(int numIn);
    void connect(int conNum, Node *nodeIn);
    Boolean isCon(int conNum);

    void printMessage(ofstream &outFile, int conNum);
    void genKCLEquation(ofstream &outFile, vector<double> &F_x, vector<double>& X, int datum, int lastnode, int nameNum, int MNAName = NA);
    void genKCLJAC(ofstream &outFile, vector<vector<double>> &JAC, vector<double>& X, int nameNum1, int nameNum2, int datum, int lastnode, int MNAName = NA);
    int genKVLEquation(ofstream &outFile, vector<double> &F_x, vector<double>& X, int datum, int lastnode);
    int genKVLJAC(ofstream &outFile, vector<vector<double>> &JAC, vector<double>& X, int datum, int lastnode);

private:
    Connectors con0, con1, con2, con3; // 器件的4个“端口”
    Component *next; // 自身预留接口，以成链
    CompType type; // 器件类型
    int compNum; // 器件编号（是在本类别中的器件编号）
    double value, temp; // 器件送值
    Model *model; // 器件对应模型
    char name[NameLength]; // 器件名称
};

class Node {// 节点
public:
    Node(int num);
    ~Node();
    int getNum();
    int getNameNum();
    int getCount();
    Node *getNext();
    Connections *getConList();
    
    void setNext(Node *nodeIn);
    void setNameNum(int numIn);
    void connect(int conNumIn, Component *compIn);

    void printMessage(ofstream &outFile);
    void genKCLEquation(ofstream &outFile, vector<double> &F_x, vector<double>& X, int datum, int lastnode, Boolean isMNA);
    void genKCLJAC(ofstream &outFile, vector<vector<double>> &JAC, vector<double>& X, int nameNum2, int datum, int lastnode, Boolean isMNA);

private:
    Node *next; // 自身预留接口，以成链
    int nodeNum, conCount; // 节点编号【这个是节点生成时的编号，并不是节点在网表中的那个编号（因为网表中的编号可能还不是按自然数顺序来的，若是的话那就差不多等价）】；节点所连组件数【用于筛选datum】
    Connections *conList; // 节点所相关的连接关系
    int nameNum; // 节点编号
};

class Model {// 模型
public:
    Model(char *nameIn, TranType typeIn, double isIn, double bfIn, double brIn, double tempIn);
    ~Model();
    Model *getNext();
    TranType getType();
    char *getName();
    double getIs();
    double getAf();
    double getAr();
    double getBf();
    double getBr();
    double getTemp();
    double getN();
    void setNext(Model *nextIn);
    
private:
    char name[NameLength]; // 模型名称
    double is, af, ar, temp, br, bf; // 模型参数
    Model *next; // 自身预留接口，以便成链
    TranType type; // 模型类型
};

class NodeHead {// 节点链表结构
public:
    NodeHead();
    ~NodeHead();
    Node *addNode();
    int getCount();
    Node *getNode(int nodeNum);

private:
    Node *nodeList; // 节点链表 
    int nodeCount; // 链上节点总数
};

class CompHead {// 器件链表结构
public:
    CompHead();
    ~CompHead();
    void addComp(Component *component);
    int getCount(CompType type);
    Component *getComp(int compNum);

private:
    Component *compList; // 器件链表
    int iCount, rCount, dCount, cCount, mCount, vSCount, iSCount, bCount; // 储存每种器件的总数目
};

class ModelHead {// 模型链表结构
public:
    ModelHead();
    ~ModelHead();
    void addModel(Model *modelIn);
    Model *getModel(char *nameIn);

private:
    Model *modelList; // 模型链表
};

class Netlist { // 网表
public:
    Netlist();
    ~Netlist();
    void setTitle(string name);
    void setAnalysisType(AnalysisType type);
    void setISIC(Boolean isIC);
    void setISNodeset(Boolean isNodeset);
    void setISOptions(Boolean isOptions);
    void insertIC(int id, double value);
    void insertNodeset(int id, double value);
    void insertOptions(string param, double value);
    void setTranStop(double stopTime);
    void setLastnode(int id);
    void setDatum(int id);
    
    ModelHead& getModelHead();
    CompHead& getCompHead();
    NodeHead& getNodeHead();
    string getTitle();
    AnalysisType getAnalysisType();
    Boolean getISIC();
    Boolean getISNodeset();
    Boolean getISOptions();
    unordered_map<int,double>& getICMap();
    unordered_map<int,double>& getNodesetMap();
    unordered_map<string,double>& getOptionsMap();
    double getTranStop();
    int getDatum();
    int getLastnode();

private:
    ModelHead modelList;
    CompHead compList;
    NodeHead nodeList;
    string title;
    AnalysisType analysisType;
    Boolean is_ic;
    Boolean is_nodeset;
    Boolean is_options;
    unordered_map<int, double> ic;
    unordered_map<int, double> nodeset;
    unordered_map<string, double> options;
    double tran_stop;
    int lastnode;
    int datum;

};

void generateMatrix(NodeHead &nodeList, CompHead &compList, ModelHead &modelList, vector<double> &F_x, vector<double> &X, vector<vector<double>> &JAC, string &outFileName, int datum, int lastnode, int step);

void parseNetList(Netlist &netlist, string &inFileName, string &outFileName);

double stripString(char *stringIn);

double calculateFe(vector<double> &X, double Is, double af, double n, int con1, int con2, int datum);

double calculateFc(vector<double> &X, double Is, double ar, double n, int con1, int con2, int datum);

double calculateFe_(vector<double> &X, double Is, double af, double n, int con1, int con2, int datum, int nameNum2);

double calculateFc_(vector<double> &X, double Is, double ar, double n, int con1, int con2, int datum, int nameNum2);

pair<int, int> getVSourceID(Connections *conList);

#endif