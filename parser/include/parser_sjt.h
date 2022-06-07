//  PARSER.H
//  Header for Parser program.

#include <ctype.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <vector>
using namespace std;

// declarations:
const double K = 1.38E-23;
const double Q = 1.60E-19;
const int NameLength = 80, BufLength = 300, NA = -1;

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

// 六种封装类
class Component;
class ComponentHead;
class Node;
class NodeHead;
class Model;
class ModelHead;

// 两种信息类
struct Connectors {// 器件“端口”类封装信息
    Flag flag; // “端口”状态
    Node *node; // 本器件该“端口”所连节点实体
    int conNum; // 本器件该“端口”所连节点编号
};
struct Connections {// 连接关系信息封装体
    Connections *next; // 下一个“连接关系”
    Component *comp; // 指向器件实体
    int conNum; // 该连接关系链上的总数？【不，这确实是“端口”编号】
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
    void genEquation(ofstream &outFile, vector<double> &F_x, vector<double>& X, int datum, int lastnode, int nameNum);
    void genJAC(ofstream &outFile, vector<vector<double>> &JAC, vector<double>& X, int nameNum1, int nameNum2, int datum, int lastnode);

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
    void genNodalEquation(ofstream &outFile, vector<double> &F_x, vector<double>& X, int datum, int lastnode);
    void genNodalJAC(ofstream &outFile, vector<vector<double>> &JAC, vector<double>& X, Node *nodePtr2, int datum, int lastnode);
    
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
    double getBf();
    double getBr();
    double getTemp();
    void setNext(Model *nextIn);
    
private:
    char name[NameLength]; // 模型名称
    double is, bf, br, temp; // 模型参数
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
    void addModel(Model *modelIn);
    Model *getModel(char *nameIn);

private:
    Model *modelList; // 模型链表
};



/*
 * 解析函数，生成节点方程矩阵F(X)与F'(X)
 * @param nodeList 节点链表
 * @param compList 器件链表
 * @param modelList 模型链表
 * @param F_x 存储节点方程的数组
 * @param X 待求解的未知向量
 * @param JAC 存储方程导数的雅克比矩阵
 * @param outFileName 输出文件名
 * @param datum 起始节点编号，默认未指定为-1
 * @param lastnode 最大节点编号
 * @param step 迭代次数
 */
void generateMatrix(NodeHead &nodeList, CompHead &compList, ModelHead &modelList, vector<double> &F_x, vector<double> &X, vector<vector<double>> &JAC, string &outFileName, int datum, int lastnode, int step);

/*
 * 解析网表文件，并将数据存储到定义的数据结构中
 * @param nodeList 节点链表，存储所有的节点
 * @param compList 器件链表，存储所有出现的器件
 * @param modelList 存储所有声明的model
 * @param datum 接地节点编号
 * @param lastnode 最大节点编号
 * @param inFileName 网表输入文件名称
 * @param outFileName 解析结果输出文件名称
 */
void parseNetList(NodeHead &nodeList, CompHead &compList, ModelHead &modelList, int &datum, int &lastnode, string &inFileName, string &outFileName);

double stripString(char *stringIn);
void printComponents(Component *compPtr);
void printNodes(Node *nodePtr, int compFlag);
char *strComponentType(Component *compPtr);