//  PARSER.CC  Ver. 1.0
//  Program to extract Nodal equations from Spice Netlist.  Ed Chan

#include "parser_v2.h"

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
void generateMatrix(NodeHead &nodeList, CompHead &compList, ModelHead &modelList, vector<double> &F_x, vector<double> &X, vector<vector<double>> &JAC, string &outFileName, int datum, int lastnode, int step) {
    vector<int> MNANodeSet1, MNANodeSet2; // 用来存储那些补偿方程的编号，比如x5-vcc是F(5)，x5-x2-v1是F(5)
    ofstream outFile;
    Node *nodePtr, *nodePtr1, *nodePtr2;
    Component *comp;

    // 打开输出文件流 【如果step不为1，那么就不去打开了】
    if(step == 1) {
        outFile.open(outFileName + "_equation.txt", ios::out);
        if(!outFile){
            cerr << "打开文件出现错误，异常退出！" << endl;
            exit(1);
        }
    }
   
    // 1. 扫描器件链表，若有电压源存在，则列KVL方程 且 记录MNA节点编号
    comp = compList.getComp(0);
    while(comp != NULL) {
        if(comp->getType() == VSource) {
            int t = comp->genKVLEquation(outFile, F_x, X, datum, lastnode);
            MNANodeSet1.push_back(t);
        }
        comp = comp->getNext();
    }

    // 2. 扫描节点列表，列各点的KCL方程
    nodePtr = nodeList.getNode(0);
    while (nodePtr != NULL) {
        if(nodePtr->getNameNum() != datum) {
            Boolean isExist = FALSE;
            for (int i = 0; i < MNANodeSet1.size();i++) {
                if(MNANodeSet1[i] == nodePtr->getNameNum())
                    isExist = TRUE;
            }
            nodePtr->genKCLEquation(outFile, F_x, X, datum, lastnode, isExist);
        }
        nodePtr = nodePtr->getNext();
    }
    outFile.close();
    outFile.clear();

    // 打开输出文件流 【如果step不为1，那么就不去打开了】
    if(step == 1) {
        outFile.open(outFileName + "_jacobian.txt", ios::out);
        if(!outFile){
            cerr << "打开文件出现错误，异常退出！" << endl;
            exit(1);
        }
    }

    // 3. 生成上面求得方程的雅可比矩阵
    // 3.1 处理含电压源的KVL方程的求导
    comp = compList.getComp(0);
    while(comp != NULL) {
        if(comp->getType() == VSource) {
            int t = comp->genKVLJAC(outFile, JAC, X, datum, lastnode);
            MNANodeSet2.push_back(t);
        }
        comp = comp->getNext();
    }
    // 3.2 处理KCL方程的求导
    nodePtr1 = nodeList.getNode(0);
    while (nodePtr1 != NULL) {
        if (nodePtr1->getNameNum() != datum) {
            nodePtr2 = nodeList.getNode(0);
            Boolean isExist = FALSE;
            for (int i = 0; i < MNANodeSet2.size();i++) { // 判断是否是已列了KVL方程的节点作为第一维度编号
                if(MNANodeSet2[i] == nodePtr1->getNameNum())
                    isExist = TRUE;
            }
            while (nodePtr2 != NULL) { // 二重循环
                if (nodePtr2->getNameNum() != datum) { 
                    nodePtr1->genKCLJAC(outFile, JAC, X, nodePtr2->getNameNum(), datum, lastnode, isExist);
                }
                nodePtr2 = nodePtr2->getNext();
            }
            if(isExist == TRUE) { // 若是的话，那么则需对其补充对自身的求导
                pair<int, int> p = getVSourceID(nodePtr1->getConList());
                if(p.second == 1) {
                    JAC[p.first + lastnode][lastnode + p.first] -= 1;
                    outFile << "JAC(" << lastnode + p.first << "," << lastnode + p.first << ") = - 1" << endl;
                } else {
                    JAC[p.first + lastnode][lastnode + p.first] += 1;
                    outFile << "JAC(" << lastnode + p.first << "," << lastnode + p.first << ") = 1" << endl;
                }     
            }
        }
        nodePtr1 = nodePtr1->getNext();
    }
    outFile.close();
    outFile.clear();

    cout << endl
         << "=============The solution of the " << step << "th iteration is completed.=============" << endl;
}

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
void parseNetList(NodeHead &nodeList, CompHead &compList, ModelHead &modelList, int& datum, int& lastnode, string& inFileName, string& outFileName) {
    ifstream inFile;
    ofstream outFile;

    // 分析过程中所需的变量
    char buf[BufLength], buf1[BufLength], buf2[BufLength], nameBuf[NameLength];
    char *bufPtr, *charPtr1, *charPtr2;
    int intBuf1, intBuf2, intBuf3, intBuf4;
    double douBuf1, douBuf2, douBuf3, douBuf4;
    CompType typeBuf;
    Component *compPtr, *compPtr1, *compPtr2;
    Node *nodePtr, *nodePtr1, *nodePtr2;
    Model *modelPtr;
    TranType TtypeBuf;

    // 2. 处理输入文件相关
    if (inFileName.empty()) {
        cerr << "Please enter the input Spice Netlist: <\"QUIT\" to exit>" << endl;
        cin >> inFileName;
        if (inFileName == "QUIT") {
            cerr << "Program Exited Abnormally!" << endl;
            exit(0);
        }
    }
    inFile.open(inFileName, ios::in);
    while (!inFile) {
        cerr << inFileName << " is an invalid input file." << endl
             << "Please enter the input Spice Netlist: <\"QUIT\" to exit>" << endl;
        cin >> inFileName;
        if (inFileName == "QUIT") {
            cerr << "Program Exited Abnormally!" << endl;
            exit(0);
        }
        inFile.open(inFileName, ios::in);
    }

    // 3. 处理输出文件相关
    if (outFileName.empty()) {
        outFileName = inFileName + ".out";
    }
    cout << "Output saved to file: " << outFileName << endl;

    // 4. 网表解析
    inFile.getline(buf, BufLength); // first line of netlist is discarded
    inFile.getline(buf, BufLength);

    // 4.1 模型遍历识别，并入链
    cout << endl
         << "==========================Models scan has started=========================" << endl;
    while (inFile.good()) {
        if ((buf == NULL) || (*buf == '\0')) {
            inFile.getline(buf, BufLength);
            continue;
        }
        strcpy(buf1, buf);
        if (!strcmp(strtok(buf1, " "), ".model")) {// strtok是一个分解字符串的函数，如果不是模型申明语句，那么跳过
            strcpy(buf2, strtok(NULL, " ")); //继续分割，赋值到buf2，此时buf2是器件的名称
            charPtr1 = strtok(NULL, " ");    //继续分割给到charPtr1，此时charPtr1是器件的类型
            if (!strcmp(charPtr1, "PNP"))    //类型处理
                TtypeBuf = PNP;
            else if (!strcmp(charPtr1, "NPN"))
                TtypeBuf = NPN;
            else if (!strcmp(charPtr1, "NMOS"))
                TtypeBuf = NMOS;
            else if (!strcmp(charPtr1, "PMOS"))
                TtypeBuf = PMOS;

            charPtr1 = strtok(NULL, " "); //继续分割，此时charPtr1后面的是晶体管的参数值
            double temp1 = NA, temp2 = NA, temp3 = NA, temp4 = NA;
            while (charPtr1 != NULL) {    //若有参数
                // 下面处理就是四种参数的值，分别对应就行。stripString是将字符串转成浮点数double的
                if ((charPtr1[0] == 'I') && (charPtr1[1] == 'S') && (charPtr1[2] == '=')) {
                    temp1 = stripString(charPtr1);
                }
                if ((charPtr1[0] == 'B') && (charPtr1[1] == 'F') && (charPtr1[2] == '=')) {
                    temp2 = stripString(charPtr1);
                }
                if ((charPtr1[0] == 'B') && (charPtr1[1] == 'R') && (charPtr1[2] == '=')) {
                    temp3 = stripString(charPtr1);
                }
                if ((charPtr1[0] == 'T') && (charPtr1[1] == 'E') && (charPtr1[4] == '=')) {
                    temp4 = stripString(charPtr1);
                }
                charPtr1 = strtok(NULL, " ");
            }
            modelPtr = new Model(buf2, TtypeBuf, temp1, temp2, temp3, temp4);
            modelList.addModel(modelPtr);
        }
        inFile.getline(buf, BufLength);
    } //模型已经扫描处理完毕
    inFile.close();
    inFile.clear();// 这点注意一下
    inFile.open(inFileName, ios::in); //现在开始扫描处理其他的信息

    // 4.2 器件遍历扫描，并入链
    cout << endl
         << "========================Components scan has started=======================" << endl;
    char model_str[9];
    inFile.getline(buf, BufLength); 
    inFile.getline(buf, BufLength);
    while (inFile.good()) {
        if ((buf == NULL) || (*buf == '\0')){
            inFile.getline(buf, BufLength);
            continue;
        }
        if (isalpha(*buf)) { //首字母是符号，因为每种组件都有对应的代号名称，如电阻是r，bjt是q，v是电压源
        //  EDIT THIS SECTION IF NEW COMPONENTS ARE ADDED!!!
        //  we could do some rearranging in this section to catch each type in order.
            switch (*buf) {
            case 'v':
            case 'V':
                typeBuf = VSource; //电压源类型
                strcpy(nameBuf, strtok(buf, " "));
                intBuf1 = atoi(strtok(NULL, " "));
                intBuf2 = atoi(strtok(NULL, " "));
                douBuf1 = atof(strtok(NULL, " "));// 该器件的属性值
                compPtr = new Component(typeBuf, douBuf1, NA, intBuf1, intBuf2, NA, NA, NULL, nameBuf);
                compList.addComp(compPtr);
                break;

            case 'i':
            case 'I':
                typeBuf = ISource; //电流源类型
                strcpy(nameBuf, strtok(buf, " "));
                intBuf1 = atoi(strtok(NULL, " "));
                intBuf2 = atoi(strtok(NULL, " "));
                douBuf1 = atof(strtok(NULL, " "));
                compPtr = new Component(typeBuf, douBuf1, NA, intBuf1, intBuf2, NA, NA, NULL, nameBuf);
                compList.addComp(compPtr);
                break;

            case 'q':
            case 'Q':
                typeBuf = BJT; // bjt类型
                strcpy(nameBuf, strtok(buf, " "));
                intBuf1 = atoi(strtok(NULL, " "));
                intBuf2 = atoi(strtok(NULL, " "));
                intBuf3 = atoi(strtok(NULL, " "));
                compPtr = new Component(typeBuf, NA, NA, intBuf1, intBuf2, intBuf3, NA,
                                        modelList.getModel(strtok(NULL, " ")), nameBuf); //找到对应的模型，给加进去
                compList.addComp(compPtr);
                break;

            case 'm':
            case 'M':
                typeBuf = MOSFET; // mos类型
                strcpy(nameBuf, strtok(buf, " "));
                intBuf1 = atoi(strtok(NULL, " "));
                intBuf2 = atoi(strtok(NULL, " "));
                intBuf3 = atoi(strtok(NULL, " "));
                intBuf4 = atoi(strtok(NULL, " "));
                compPtr = new Component(typeBuf, NA, NA, intBuf1, intBuf2, intBuf3, intBuf4,
                                        modelList.getModel(strtok(NULL, " ")), nameBuf);
                compList.addComp(compPtr);
                break;

            case 'r':
            case 'R':
                typeBuf = Resistor; //电阻类型
                strcpy(nameBuf, strtok(buf, " "));
                intBuf1 = atoi(strtok(NULL, " "));
                intBuf2 = atoi(strtok(NULL, " "));
                douBuf1 = atof(strtok(NULL, " "));
                compPtr = new Component(typeBuf, douBuf1, NA, intBuf1, intBuf2, NA, NA, NULL, nameBuf);
                compList.addComp(compPtr);
                break;

            case 'd':
            case 'D':
                typeBuf = Diode; //就是只有这种类型的组件有tempin送值
                strcpy(nameBuf, strtok(buf, " "));
                intBuf1 = atoi(strtok(NULL, " "));
                intBuf2 = atoi(strtok(NULL, " "));
                charPtr1 = strtok(NULL, " ");
                while (charPtr1 != NULL) {
                    if ((charPtr1[0] == 'I') && (charPtr1[1] == 'S') && (charPtr1[2] == '=')) {
                        douBuf1 = stripString(charPtr1);
                    }
                    if ((charPtr1[0] == 'T') && (charPtr1[1] == 'E') && (charPtr1[4] == '=')) {
                        douBuf2 = stripString(charPtr1);
                    }
                    charPtr1 = strtok(NULL, " ");
                }
                compPtr = new Component(typeBuf, douBuf1, douBuf2, intBuf1, intBuf2, NA, NA, NULL, nameBuf);
                compList.addComp(compPtr);
                break;

            case 'c':
            case 'C':
                typeBuf = Capacitor;
                strcpy(nameBuf, strtok(buf, " "));
                intBuf1 = atoi(strtok(NULL, " "));
                intBuf2 = atoi(strtok(NULL, " "));
                douBuf1 = atof(strtok(NULL, " "));
                compPtr = new Component(typeBuf, douBuf1, NA, intBuf1, intBuf2, NA, NA, NULL, nameBuf);
                compList.addComp(compPtr);
                break;

            case 'l':
            case 'L':
                typeBuf = Inductor;
                strcpy(nameBuf, strtok(buf, " "));
                intBuf1 = atoi(strtok(NULL, " "));
                intBuf2 = atoi(strtok(NULL, " "));
                douBuf1 = atof(strtok(NULL, " "));
                compPtr = new Component(typeBuf, douBuf1, NA, intBuf1, intBuf2, NA, NA, NULL, nameBuf);
                compList.addComp(compPtr);
                break;
            }
        }
        inFile.getline(buf, BufLength);
    }
    inFile.close();
    inFile.clear();

    cout << endl
         << "=============The contents of the netlist have been scanned================" << endl
         << endl;

    /************************************************************/
    /** 对于起始节点编号的确定是一个问题。常见的是接地点为0号节点， **/
    /** 也即为第一个节点。但似乎有约定是所连器件数最多的节点作为公  **/
    /** 共节点，即接地节点。                                    **/
    /************************************************************/

    // 4.3 依据器件链表，进行节点遍历扫描，节点成链并与器件建立联系
    cout << endl
         << "=================The node generation process has started==================" << endl;
    compPtr1 = compList.getComp(0);
    while (compPtr1 != NULL) {
        for (int b = 0; b < 3; b++) { // TODO：遇到4个端口的器件需注意这里的改变
            // 验证端口b是否被遍历过了 &&  端口所连节点编号存在
            if ((!compPtr1->isCon(b)) && (compPtr1->getConVal(b) != NA)) {
                intBuf1 = compPtr1->getConVal(b);
                //【这边应该有一个去重的操作，不然会重复创建节点。这里就不用去重了，因为下面人家也是直接扫描了，组件connect的时候会置SET】
                nodePtr1 = nodeList.addNode();  //新生成节点，并在每个节点中记录当前已有多少个节点【其实是节点编号】，同样在头结点链表中也有记录
                nodePtr1->setNameNum(intBuf1);  // ~> naming the node as in the netlist file
                compPtr1->connect(b, nodePtr1); // ~> connecting the 'connector' of component to the node
                nodePtr1->connect(b, compPtr1); // ~> connecting the 'connection' of the node to the component

                // now search and connect all other appropriate connectors to this node.
                // error checking should be added to prevent duplicated, or skipped connectors.
                // compPtr2 = compPtr1->getNext();
                compPtr2 = compPtr1->getNext();
                while (compPtr2 != NULL) {
                    for (int c = 0; c < 3; c++) {
                        if (compPtr2->getConVal(c) == intBuf1) { 
                            compPtr2->connect(c, nodePtr1);
                            nodePtr1->connect(c, compPtr2);
                            break;
                        }
                    }
                    compPtr2 = compPtr2->getNext();
                }
            }
        }
        compPtr1 = compPtr1->getNext();
    }

    // 5. 确认传入的起始节点编号可用、正确；如果没有传入起始节点编号，那么找到连接器件数最多的节点作为接地节点
    if (datum != NA) {
        Boolean check = FALSE;
        nodePtr = nodeList.getNode(0);
        while (nodePtr != NULL) {
            if (nodePtr->getNameNum() == datum)
                check = TRUE;
            nodePtr = nodePtr->getNext();
        }
        if (check == FALSE) {
            cerr << "Datum value invalid!" << endl
                 << "PROGRAM EXITED ABNORMALLY!" << endl;
            exit(0);
        }
    } else {
        nodePtr = nodeList.getNode(0);
        nodePtr1 = nodePtr->getNext();
        while (nodePtr1 != NULL) {
            if (nodePtr1->getCount() > nodePtr->getCount())
                nodePtr = nodePtr1;
            nodePtr1 = nodePtr1->getNext();
        }
        datum = nodePtr->getNameNum(); // datum是Count数最大的那个
    }

    // 6. 找到节点编号中最大的
    // 注意：不能简单地通过nodeList中的nodeCount就得出编号，因为对于网表描述的不一，有的0接地，有的不存在0
    nodePtr = nodeList.getNode(0); //~> getting the pointer to the first node, pointed by 'headNode'
    lastnode = nodePtr->getNameNum();
    while (nodePtr != NULL) {
        lastnode = (nodePtr->getNameNum() > lastnode) ? nodePtr->getNameNum() : lastnode;
        nodePtr = nodePtr->getNext();
    }

    // 7. 输出解析内容，检查是否正确解析
    outFile.open(outFileName + "_parser.txt", ios::out);
    if(!outFile.is_open()) {
        cerr << "Failed to open " << outFileName + "_parser.txt" << '\n';
        exit(0);
    }
    outFile << "datum = " << datum << "        lastnode = " << lastnode << endl;
    nodePtr = nodeList.getNode(0);
    while(nodePtr != NULL){
        outFile << "节点"<<nodePtr->getNameNum()<<"        所连器件数为："<<nodePtr->getCount() << endl;
        nodePtr->printMessage(outFile);
        nodePtr = nodePtr->getNext();
    }
    outFile.close();
    outFile.clear();

    cout << endl
         << "====================The netlist is parsed completely======================" << endl;

    return;
}

// 将传入的字符串解析为double型数字
double stripString(char *stringIn) {
    char buf[BufLength], buf2[BufLength];
    int a, b;
    strcpy(buf, stringIn);
    for (a = 0; buf[a] != '='; a++);
    a++;
    for (b = 0; buf[a] != '\0'; b++, a++)
        buf2[b] = buf[a];
    buf2[b] = '\0';
    return atof(buf2); //转成浮点数
}

// 计算bjt模型的fe函数值
double calculateFe(vector<double> &X, double Is, double af, double n, int con2, int con1, int datum) {
    if (con2 != datum && con1 != datum)
        return (-Is) / af * (exp(-n * (X[con2] - X[con1])) - 1);
    else if (con2 != datum && con1 == datum)
        return (-Is) / af * (exp(-n * X[con2]) - 1);
    else if (con2 == datum && con1 != datum)
        return (-Is) / af * (exp(n * X[con1]) - 1);
    else
        return 0;
}

// 计算bjt模型的fc函数值
double calculateFc(vector<double> &X, double Is, double ar, double n, int con0, int con1, int datum) {
    if (con0 != datum && con1 != datum)
        return (-Is) / ar * (exp(-n * (X[con0] - X[con1])) - 1);
    else if (con0 != datum && con1 == datum)
        return (-Is) / ar * (exp(-n * X[con0]) - 1);
    else if (con0 == datum && con1 != datum)
        return (-Is) / ar * (exp(n * X[con1]) - 1);
    else
        return 0;
}

// 计算bjt模型的fe函数的导数值
double calculateFe_(vector<double> &X, double Is, double af, double n, int con2, int con1, int datum, int nameNum2) {
    double temp = 0;
    if (con2 != datum && con1 != datum) {
        temp = (-Is) * (-n) / af * exp(-n * (X[con2] - X[con1]));
        temp = nameNum2 == con1 ? (-temp) : temp;
    }
    else if (con2 != datum && con1 == datum)
        temp = (-Is) * (-n) / af * exp(-n * X[con2]);
    else if (con2 == datum && con1 != datum)
        temp = (-Is) * n / af * exp(n * X[con1]);
    return temp;
}

// 计算bjt模型的fc函数的导数值
double calculateFc_(vector<double> &X, double Is, double ar, double n, int con0, int con1, int datum, int nameNum2) {
    double temp = 0;
    if (con0 != datum && con1 != datum) {
        temp = (-Is) * (-n) / ar * exp(-n * (X[con0] - X[con1]));
        temp = nameNum2 == con1 ? (-temp) : temp;
    }
    else if (con0 != datum && con1 == datum)
        temp = (-Is) * (-n) / ar * exp(-n * X[con0]);    
    else if (con0 == datum && con1 != datum)
        temp = (-Is) * n / ar * exp(n * X[con1]);
    return temp;
}

/*
 * 获得某节点的连接关系中，电压源名字与所连端口编号
 * @param conList 传入的连接关系
 */
pair<int, int> getVSourceID(Connections *conList) {
    while(conList != NULL) {
        if(conList->comp->getType() == VSource)
            return {conList->comp->getcompNum(), conList->conNum};
        conList = conList->next;
    }
    return {NA, NA};
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                     //
//               上面是“启动”函数，下面是本parser中采用的数据结构相关的类中成员函数的实现                    //
//                                                                                                     //
/////////////////////////////////////////////////////////////////////////////////////////////////////////


Component::Component(CompType typeIn, double valueIn = NA, double tempIn = NA,
                     int con0In = NA, int con1In = NA, int con2In = NA, int con3In = NA,
                     Model *modelIn = NULL, char *nameIn = NULL) {

    type = typeIn;
    con0.conNum = con0In;
    con1.conNum = con1In;
    con2.conNum = con2In;
    con3.conNum = con3In;
    con0.flag = UNSET;
    con1.flag = UNSET;
    con2.flag = UNSET;
    con3.flag = UNSET;
    value = valueIn;
    temp = tempIn;
    next = NULL;
    model = modelIn;
    strcpy(name, nameIn);
}

Component::~Component(){};

/*
 * 将节点实体连接到器件的对应端口上
 * @param conNum 要连接到器件的“端口”号
 * @param nodeIn 要连接到器件上的节点实体
 */
void Component::connect(int conNum, Node *nodeIn) {
    if (conNum == 0) {
        con0.node = nodeIn;
        con0.flag = SET;
    }
    if (conNum == 1) {
        con1.node = nodeIn;
        con1.flag = SET;
    }
    if (conNum == 2) {
        con2.node = nodeIn;
        con2.flag = SET;
    }
    if (conNum == 3) {
        con3.node = nodeIn;
        con3.flag = SET;
    }
}

// 获取器件类型
CompType Component::getType() {
    return type;
}

int Component::getNum() {
    return compNum;
}

// 获取器件编号
int Component::getcompNum() {
    return compNum;
}

// 获取器件链表中的下一项
Component *Component::getNext() {
    return next;
}

// 获得器件的值
double Component::getVal() {
    return value;
}

// 添加组件到链上，并完成链接关系
void Component::setNext(Component *nextIn) {
    next = nextIn;
}

// 设置器件编号
void Component::setNum(int numIn) {
    compNum = numIn;
}

// 获得本器件某端口所连节点编号
int Component::getConVal(int conNum) {
    int rtVal;
    if (conNum == 0)
        rtVal = con0.conNum;
    if (conNum == 1)
        rtVal = con1.conNum;
    if (conNum == 2)
        rtVal = con2.conNum;
    if (conNum == 3)
        rtVal = con3.conNum;
    return rtVal;
}

// 判断端口是否可用
Boolean Component::isCon(int conNum) {
    Boolean rtVal;
    if (conNum == 0)
        rtVal = (con0.flag == SET) ? TRUE : FALSE;
    if (conNum == 1)
        rtVal = (con1.flag == SET) ? TRUE : FALSE;
    if (conNum == 2)
        rtVal = (con2.flag == SET) ? TRUE : FALSE;
    if (conNum == 3)
        rtVal = (con3.flag == SET) ? TRUE : FALSE;
    return rtVal;
}

// 从端口获取节点实体
Node *Component::getNode(int conNum) {
    switch (conNum){
    case 0:
        return con0.node;
    case 1:
        return con1.node;
    case 2:
        return con2.node;
    case 3:
        return con3.node;
    }
    return NULL;
}

// 从端口获取节点数目【这个没啥用，不用关注】
int Component::getNodeNum(int conNum) {
    switch (conNum) {
    case 0:
        return con0.node->getNum();
    case 1:
        return con1.node->getNum();
    case 2:
        return con2.node->getNum();
    case 3:
        return con3.node->getNum();
    }
    return NA;
}

// 获取器件名称
char *Component::getName() {
    return name;
}

// 打印器件的信息
void Component::printMessage(ofstream &outFile, int conNum){
    switch (type) {
    case BJT:
        outFile << "        编号：" << getcompNum() << "    类型："
                << "BJT"
                << " 连接端口：" << conNum << "    名称：" << name << endl;
        outFile << "        value: IS = " << model->getIs() << "  AF = " << model->getAf() << "  AR = " << model->getAr() << "  N = " << model->getN() << endl;
        break;
    case VSource:
        outFile << "        编号：" << getcompNum() << "    类型："
                << "VSource"
                << " 连接端口：" << conNum << "    名称：" << name << endl;
        outFile << "        value: " << value << endl;
        break;
    case Resistor:
        outFile << "        编号：" << getcompNum() << "    类型："
                << "Resistor"
                << " 连接端口：" << conNum << "    名称：" << name << endl;
        outFile << "        value: " << value << endl;
        break;
    case ISource:
        outFile << "        编号：" << getcompNum() << "    类型："
                << "ISource"
                << " 连接端口：" << conNum << "    名称：" << name << endl;
        outFile << "        value: " << value << endl;    
        break;
    }
}

/*******************************************************************************************************************************************************************/

/*
 * 生成本组件相关的方程（KCL）
 * @param outFile 输出文件流
 * @param F_x F(X)方程组
 * @param X 待求解向量
 * @param datum 接地节点编号
 * @param lastnode 最大节点编号
 * @param nameNum 节点编号（为该节点生成方程）
 * @param MNAName 记录未知电流变量编号
 */
void Component::genKCLEquation(ofstream &outFile, vector<double> &F_x, vector<double>& X, int datum, int lastnode, int nameNum, int MNAName) {
    int actualName = MNAName == NA ? nameNum : MNAName;
    switch (type) {
    case MOSFET:
        //这里存疑？是否MOSFET与BJT是一个特性的东西
    case BJT: {
        double is = model->getIs(), ar = model->getAr(), af = model->getAf(), n = model->getN();
        double fc = calculateFc(X, is, ar, n, con0.node->getNameNum(), con1.node->getNameNum(), datum);
        double fe = calculateFe(X, is, af, n, con2.node->getNameNum(), con1.node->getNameNum(), datum);
        if ((con0.node->getNameNum() == nameNum) && (model->getType() == NPN)) {
            F_x[actualName] += fc - af * fe;
            outFile << " + " << name << "_Ic";
        }
        else if ((con2.node->getNameNum() == nameNum) && (model->getType() == NPN)) {
            F_x[actualName] += fe - ar * fc;
            outFile << " + " << name << "_Ie";     
        }
        else if((con1.node->getNameNum() == nameNum) && (model->getType() == NPN)) {
            F_x[actualName] -= fc - af * fe + fe - ar * fc;
            outFile << " - " << name << "_Ie"
                    << " - " << name << "_Ic";
        }
    }
        // TODO: 其他情况是知识完善后，再做补充处理
        break;

    case VSource: // KCL
        if(con0.node->getNameNum() == nameNum) {
            F_x[actualName] += X[lastnode + compNum];
            outFile << " + x(" << lastnode + compNum << ")";
        }
        else if(con1.node->getNameNum() == nameNum) {
            F_x[actualName] -= X[lastnode + compNum];
            outFile << " - x(" << lastnode + compNum << ")";
        }
        break;

    case ISource:
        if (con0.node->getNameNum() == nameNum) {
            F_x[actualName] += value;
            outFile << " + " << name;
        } else if (con1.node->getNameNum() == nameNum) {
            F_x[actualName] -= value;
            outFile << " - " << name;
        }
        break;

    case Diode:
        // TODO：待完成该情况下处理
        break;

    case Resistor:
        if (con0.node->getNameNum() == nameNum) {
            // F_x[actualName] += (X[nameNum] - X[con1.node->getNameNum()]) / value;
            outFile << " + (";
            if (con0.node->getNameNum() != datum) {
                F_x[actualName] += X[con0.node->getNameNum()] / value;
                outFile << " x(" << con0.node->getNameNum() << ')';
            }
            if (con1.node->getNameNum() != datum) {
                F_x[actualName] -= X[con1.node->getNameNum()] / value;
                outFile << " - x(" << con1.node->getNameNum() << ')';
            }  
            outFile << " ) / " << name;
        }
        else if (con1.node->getNameNum() == nameNum) {
            // F_x[actualName] += (X[nameNum] - X[con0.node->getNameNum()]) / value;
            outFile << " + (";
            if (con1.node->getNameNum() != datum) {
                F_x[actualName] += X[con1.node->getNameNum()] / value;
                outFile << " x(" << con1.node->getNameNum() << ')';
            }
            if (con0.node->getNameNum() != datum) {
                F_x[actualName] -= X[con0.node->getNameNum()] / value;
                outFile << " - x(" << con0.node->getNameNum() << ')';
            }   
            outFile << " ) / " << name;
        }
        break;

    case Capacitor:
        // TODO：待完成该情况下处理
        break;

    case Inductor:
        // TODO：待完成该情况下处理
        break;
    }
    return;
}

/*
 * 生成本组件相关的方程（KVL），并返回“方程”的编号
 * @param outFile 输出文件流
 * @param F_x F(X)方程组
 * @param X 待求解向量
 * @param datum 接地节点编号
 * @param lastnode 最大节点编号
 * @param nameNum 节点编号（为该节点生成方程）
 * @param MNAName 记录未知电流变量编号
 */
int Component::genKVLEquation(ofstream &outFile, vector<double> &F_x, vector<double>& X, int datum, int lastnode) {
    if(con0.conNum != datum && con1.conNum != datum) {
        F_x[con0.conNum] += X[con0.conNum] - X[con1.conNum] - value;
        outFile << "F(" << con0.conNum << ") = x(" << con0.conNum << ") - x(" << con1.conNum << ") - " << name << endl;
        return con0.conNum;
    }
    else if(con0.conNum != datum && con1.conNum == datum) {
        F_x[con0.conNum] += X[con0.conNum] - value;
        outFile << "F(" << con0.conNum << ") = x(" << con0.conNum << ") - " << name << endl;
        return con0.conNum;
    }
    else if(con0.conNum == datum && con1.conNum != datum) {
        F_x[con1.conNum] -= X[con1.conNum] + value;
        outFile << "F(" << con1.conNum << ") = - x(" << con1.conNum << ") - " << name << endl;
        return con1.conNum;
    }
    return NA;
}

/*
 * 生成本组件相关的方程（KVL）的导数，并返回系列“方程”（对应JAC的一维）编号
 * @param outFile 输出文件流
 * @param JAC F(X)方程组的雅可比矩阵
 * @param X 待求解向量
 * @param datum 接地节点编号
 * @param lastnode 最大节点编号
 */
int Component::genKVLJAC(ofstream &outFile, vector<vector<double>> &JAC, vector<double>& X, int datum, int lastnode) {
    if(con0.conNum != datum && con1.conNum != datum) {
        JAC[con0.conNum][con0.conNum] += 1;
        JAC[con0.conNum][con1.conNum] -= 1;
        outFile << "JAC(" << con0.conNum << "," << con0.conNum << ") = 1" << endl;
        outFile << "JAC(" << con0.conNum << "," << con1.conNum << ") = - 1" << endl;
        return con0.conNum;
    }
    else if(con0.conNum != datum && con1.conNum == datum) {
        JAC[con0.conNum][con0.conNum] += 1;
        outFile << "JAC(" << con0.conNum << "," << con0.conNum << ") = 1" << endl;
        return con0.conNum;
    }
    else if(con0.conNum == datum && con1.conNum != datum) {
        JAC[con1.conNum][con1.conNum] -= 1;
        outFile << "JAC(" << con1.conNum << "," << con1.conNum << ") = - 1" << endl;
        return con1.conNum;
    }
    return NA;
}

/*
 * 生成本组件相关的方程（KCL）的导数
 * @param outFile 输出文件流
 * @param JAC F(X)方程组的雅可比矩阵
 * @param X 待求解向量
 * @param nameNum1 一维编号
 * @param nameNum2 二维编号
 * @param datum 接地节点编号
 * @param lastnode 最大节点编号
 * @param MNAName 记录未知电流变量编号
 */
void Component::genKCLJAC(ofstream &outFile, vector<vector<double>> &JAC, vector<double>& X, int nameNum1, int nameNum2, int datum, int lastnode, int MNAName) {
    int actualName = MNAName == NA ? nameNum1 : MNAName;
    switch (type) {
    case MOSFET:
        std::cout << "genJAC for MOSFETs not implemented" << std::endl
                  << "PROGRAM ENDED ABNORMALLY!" << std::endl;
        exit(0);
        break;

    case BJT: {
        double is = model->getIs(), ar = model->getAr(), af = model->getAf(), n = model->getN();
        if ((con0.node->getNameNum() == nameNum1) && (model->getType() == NPN) && (con0.node->getNameNum() == nameNum2)) {
            // JAC[actualName][nameNum2] += (-1.0 * 1e-16) * (-38.78) / 0.5 * exp(38.78 * (X[con0.node->getNameNum()] - X[con1.node->getNameNum()]));
            outFile << " + IS * N / ar * exp( - N * ( ";
            if (con0.node->getNameNum() != datum)
                outFile << "x(" << con0.node->getNameNum() << ")";
            if (con1.node->getNameNum() != datum)
                outFile << " - x(" << con1.node->getNameNum() << ")";
            outFile << " ) )";

            // 这里也像下面一样，注释掉一些东西，不过删除了

            JAC[actualName][nameNum2] += calculateFc_(X, is, ar, n, con0.node->getNameNum(), con1.node->getNameNum(), datum, nameNum2);
        }
        else if ((con0.node->getNameNum() == nameNum1) && (model->getType() == NPN) && (con2.node->getNameNum() == nameNum2)) {
            // JAC[actualName][nameNum2] += (1.0 * 1e-16) * (-38.78) * exp(38.78 * (X[con2.node->getNameNum()] - X[con1.node->getNameNum()]));
            outFile << " - IS * N * exp( - N * ( ";
            if (con2.node->getNameNum() != datum)
                outFile << "x(" << con2.node->getNameNum() << ")";
            if (con1.node->getNameNum() != datum)
                outFile << " - x(" << con1.node->getNameNum() << ")";
            outFile << " ) )";
            // if (con2.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
            //     JAC[actualName][nameNum2] += (1.0 * 1e-16) * (-38.78) * exp(38.78 * (X[con2.node->getNameNum()] - X[con1.node->getNameNum()]));
            // } else if (con2.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
            //     JAC[actualName][nameNum2] += (1.0 * 1e-16) * (-38.78) * exp(38.78 * X[con2.node->getNameNum()]);
            // } else {
            //     JAC[actualName][nameNum2] += (1.0 * 1e-16) * (-38.78) * exp(38.78 * (-X[con1.node->getNameNum()]));
            // }
            JAC[actualName][nameNum2] -= af * calculateFe_(X, is, af, n, con2.node->getNameNum(), con1.node->getNameNum(), datum, nameNum2);
        }
        else if ((con0.node->getNameNum() == nameNum1) && (model->getType() == NPN) && (con1.node->getNameNum() == nameNum2)) {
            // double temp1 = (-1.0 * 1e-16) * (-38.78) * exp(38.78 * (X[con2.node->getNameNum()] - X[con1.node->getNameNum()]));
            // double temp2 = (-1.0 * 1e-16) * (-38.78) / 0.5 * exp(38.78 * (X[con0.node->getNameNum()] - X[con1.node->getNameNum()]));
            // JAC[actualName][nameNum2] += temp1 - temp2;
            outFile << " + IS * N * exp( - N * ( ";
            if (con2.node->getNameNum() != datum)
                outFile << "x(" << con2.node->getNameNum() << ")";
            if (con1.node->getNameNum() != datum)
                outFile << " - x(" << con1.node->getNameNum() << ")";
            // if (con2.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
            //     JAC[actualName][nameNum2] += (-1.0 * 1e-16) * (-38.78) * exp(38.78 * (X[con2.node->getNameNum()] - X[con1.node->getNameNum()]));
            // } else if (con2.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
            //     JAC[actualName][nameNum2] += (-1.0 * 1e-16) * (-38.78) * exp(38.78 * X[con2.node->getNameNum()]);
            // } else {
            //     JAC[actualName][nameNum2] += (-1.0 * 1e-16) * (-38.78) * exp(38.78 * (-X[con1.node->getNameNum()]));
            // }
            JAC[actualName][nameNum2] -= af * calculateFe_(X, is, af, n, con2.node->getNameNum(), con1.node->getNameNum(), datum, nameNum2);
            outFile << " ) )";
            outFile << " - ";
            outFile << "IS * N / ar * exp( - N * ( ";
            if (con0.node->getNameNum() != datum)
                outFile << "x(" << con0.node->getNameNum() << ")";
            if (con1.node->getNameNum() != datum)
                outFile << " - x(" << con1.node->getNameNum() << ")";
            outFile << " ) )";
            // if (con0.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
            //     JAC[actualName][nameNum2] -= (-1.0 * 1e-16) * (-38.78) / 0.5 * exp(38.78 * (X[con0.node->getNameNum()] - X[con1.node->getNameNum()]));
            // } else if (con0.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
            //     JAC[actualName][nameNum2] -= (-1.0 * 1e-16) * (-38.78) / 0.5 * exp(38.78 * X[con0.node->getNameNum()]);
            // } else {
            //     JAC[actualName][nameNum2] -= (-1.0 * 1e-16) * (-38.78) / 0.5 * exp(38.78 * (-X[con1.node->getNameNum()]));
            // }
            JAC[actualName][nameNum2] += calculateFc_(X, is, ar, n, con0.node->getNameNum(), con1.node->getNameNum(), datum, nameNum2);
        }
        else if ((con2.node->getNameNum() == nameNum1) && (model->getType() == NPN) && (con0.node->getNameNum() == nameNum2)) {
            // JAC[actualName][nameNum2] -= (-1.0 * 1e-16) * (-38.78) * exp(38.78 * (X[con0.node->getNameNum()] - X[con1.node->getNameNum()]));
            outFile << " - IS * N * exp( - N * ( ";
            if (con0.node->getNameNum() != datum)
                outFile << "x(" << con0.node->getNameNum() << ")";
            if (con1.node->getNameNum() != datum)
                outFile << " - x(" << con1.node->getNameNum() << ")";
            outFile << " ) )";
            // if (con0.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
            //     JAC[actualName][nameNum2] -= (-1.0 * 1e-16) * (-38.78) * exp(38.78 * (X[con0.node->getNameNum()] - X[con1.node->getNameNum()]));
            // } else if (con0.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
            //     JAC[actualName][nameNum2] -= (-1.0 * 1e-16) * (-38.78) * exp(38.78 * X[con0.node->getNameNum()]);
            // } else {
            //     JAC[actualName][nameNum2] -= (-1.0 * 1e-16) * (-38.78) * exp(38.78 * (-X[con1.node->getNameNum()]));
            // }
            JAC[actualName][nameNum2] -= ar * calculateFc_(X, is, ar, n, con0.node->getNameNum(), con1.node->getNameNum(), datum, nameNum2);
        }
        else if ((con2.node->getNameNum() == nameNum1) && (model->getType() == NPN) && (con2.node->getNameNum() == nameNum2)) {
            // JAC[actualName][nameNum2] += (-1.0 * 1e-16) * (-38.78) / 0.99 * exp(38.78 * (X[con2.node->getNameNum()] - X[con1.node->getNameNum()]));
            outFile << " + IS * N / af * exp( - N * ( ";
            if (con2.node->getNameNum() != datum)
                outFile << "x(" << con2.node->getNameNum() << ')';
            if (con1.node->getNameNum() != datum)
                outFile << " - x(" << con1.node->getNameNum() << ')';
            outFile << " ) )";
            // if (con2.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
            //     JAC[actualName][nameNum2] += (-1.0 * 1e-16) * (-38.78) / 0.99 * exp(38.78 * (X[con2.node->getNameNum()] - X[con1.node->getNameNum()]));
            // } else if (con2.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
            //     JAC[actualName][nameNum2] += (-1.0 * 1e-16) * (-38.78) / 0.99 * exp(38.78 * X[con2.node->getNameNum()]);
            // } else {
            //     JAC[actualName][nameNum2] += (-1.0 * 1e-16) * (-38.78) / 0.99 * exp(38.78 * (-X[con1.node->getNameNum()]));
            // }
            JAC[actualName][nameNum2] += calculateFe_(X, is, af, n, con2.node->getNameNum(), con1.node->getNameNum(), datum, nameNum2);
        }
        else if ((con2.node->getNameNum() == nameNum1) && (model->getType() == NPN) && (con1.node->getNameNum() == nameNum2)) {
            // double temp2 = (-1.0 * 1e-16) * (-38.78) * exp(38.78 * (X[con0.node->getNameNum()] - X[con1.node->getNameNum()]));
            // double temp1 = (-1.0 * 1e-16) * (-38.78) / 0.99 * exp(38.78 * (X[con2.node->getNameNum()] - X[con1.node->getNameNum()]));
            // JAC[actualName][nameNum2] += temp2 - temp1;
            outFile << " - IS * N / af * exp( - N * ( ";
            if (con2.node->getNameNum() != datum)
                outFile << "x(" << con2.node->getNameNum() << ")";
            if (con1.node->getNameNum() != datum)
                outFile << " - x(" << con1.node->getNameNum() << ")";
            // if (con2.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
            //     JAC[actualName][nameNum2] -= (-1.0 * 1e-16) * (-38.78) / 0.99 * exp(38.78 * (X[con2.node->getNameNum()] - X[con1.node->getNameNum()]));
            // } else if (con2.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
            //     JAC[actualName][nameNum2] -= (-1.0 * 1e-16) * (-38.78) / 0.99 * exp(38.78 * X[con2.node->getNameNum()]);
            // } else {
            //     JAC[actualName][nameNum2] -= (-1.0 * 1e-16) * (-38.78) / 0.99 * exp(38.78 * (-X[con1.node->getNameNum()]));
            // }
            JAC[actualName][nameNum2] += calculateFe_(X, is, af, n, con2.node->getNameNum(), con1.node->getNameNum(), datum, nameNum2);
            outFile << " ) )";
            outFile << " + ";
            outFile << "IS * N * exp( - N * ( ";
            if (con0.node->getNameNum() != datum)
                outFile << "x(" << con0.node->getNameNum() << ")";
            if (con1.node->getNameNum() != datum)
                outFile << " - x(" << con1.node->getNameNum() << ")";
            outFile << " ) )";
            // if (con0.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
            //     JAC[actualName][nameNum2] += (-1.0 * 1e-16) * (-38.78) * exp(38.78 * (X[con0.node->getNameNum()] - X[con1.node->getNameNum()]));
            // } else if (con0.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
            //     JAC[actualName][nameNum2] += (-1.0 * 1e-16) * (-38.78) * exp(38.78 * X[con0.node->getNameNum()]);
            // } else {
            //     JAC[actualName][nameNum2] += (-1.0 * 1e-16) * (-38.78) * exp(38.78 * (-X[con1.node->getNameNum()]));
            // }
            JAC[actualName][nameNum2] -= ar * calculateFc_(X, is, ar, n, con0.node->getNameNum(), con1.node->getNameNum(), datum, nameNum2);
        }
        else if ((con1.node->getNameNum() == nameNum1) && (model->getType() == NPN) && (con0.node->getNameNum() == nameNum2)) {
            // double temp1 = (-1.0 * 1e-16) * (-38.78) * exp(38.78 * (X[con0.node->getNameNum()] - X[con1.node->getNameNum()]));
            // double temp2 = (-1.0 * 1e-16) * (-38.78) / 0.5 * exp(38.78 * (X[con0.node->getNameNum()] - X[con1.node->getNameNum()]));
            // JAC[actualName][nameNum2] += temp1 - temp2;
            outFile << " + IS * N * exp( - N * ( ";
            if (con0.node->getNameNum() != datum)
                outFile << "x(" << con0.node->getNameNum() << ")";
            if (con1.node->getNameNum() != datum)
                outFile << " - x(" << con1.node->getNameNum() << ")";
            // if (con0.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
            //     JAC[actualName][nameNum2] += (-1.0 * 1e-16) * (-38.78) * exp(38.78 * (X[con0.node->getNameNum()] - X[con1.node->getNameNum()]));
            // } else if (con0.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
            //     JAC[actualName][nameNum2] += (-1.0 * 1e-16) * (-38.78) * exp(38.78 * X[con0.node->getNameNum()]);
            // } else {
            //     JAC[actualName][nameNum2] += (-1.0 * 1e-16) * (-38.78) * exp(38.78 * (-X[con1.node->getNameNum()]));
            // }
            JAC[actualName][nameNum2] += ar * calculateFc_(X, is, ar, n, con0.node->getNameNum(), con1.node->getNameNum(), datum, nameNum2);
            outFile << " ) )";
            outFile << " - ";
            outFile << "IS * N / ar * exp( - N * ( ";
            if (con0.node->getNameNum() != datum)
                outFile << "x(" << con0.node->getNameNum() << ")";
            if (con1.node->getNameNum() != datum)
                outFile << " - x(" << con1.node->getNameNum() << ")";
            outFile << " ) ) ";
            // if (con0.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
            //     JAC[actualName][nameNum2] -= (-1.0 * 1e-16) * (-38.78) / 0.5 * exp(38.78 * (X[con0.node->getNameNum()] - X[con1.node->getNameNum()]));
            // } else if (con0.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
            //     JAC[actualName][nameNum2] -= (-1.0 * 1e-16) * (-38.78) / 0.5 * exp(38.78 * X[con0.node->getNameNum()]);
            // } else {
            //     JAC[actualName][nameNum2] -= (-1.0 * 1e-16) * (-38.78) / 0.5 * exp(38.78 * (-X[con1.node->getNameNum()]));
            // }
            JAC[actualName][nameNum2] -= calculateFc_(X, is, ar, n, con0.node->getNameNum(), con1.node->getNameNum(), datum, nameNum2);
        }
        else if ((con1.node->getNameNum() == nameNum1) && (model->getType() == NPN) && (con2.node->getNameNum() == nameNum2)) {
            // double temp1 = (-1.0 * 1e-16) * (-38.78) / 0.99 * exp(38.78 * (X[con2.node->getNameNum()] - X[con1.node->getNameNum()]));
            // double temp2 = (-1.0 * 1e-16) * (-38.78) * exp(38.78 * (X[con2.node->getNameNum()] - X[con1.node->getNameNum()]));
            // JAC[actualName][nameNum2] += temp2 - temp1;
            outFile << " + IS * N * exp( - N * ( ";
            if (con2.node->getNameNum() != datum)
                outFile << "x(" << con2.node->getNameNum() << ")";
            if (con1.node->getNameNum() != datum)
                outFile << " - x(" << con1.node->getNameNum() << ")";
            // if (con2.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
            //     JAC[actualName][nameNum2] += (-1.0 * 1e-16) * (-38.78) * exp(38.78 * (X[con2.node->getNameNum()] - X[con1.node->getNameNum()]));
            // } else if (con2.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
            //     JAC[actualName][nameNum2] += (-1.0 * 1e-16) * (-38.78) * exp(38.78 * X[con2.node->getNameNum()]);
            // } else {
            //     JAC[actualName][nameNum2] += (-1.0 * 1e-16) * (-38.78) * exp(38.78 * (-X[con1.node->getNameNum()]));
            // }
            JAC[actualName][nameNum2] += af * calculateFe_(X, is, af, n, con2.node->getNameNum(), con1.node->getNameNum(), datum, nameNum2);
            outFile << " ) )";
            outFile << " - ";
            outFile << "IS * N / af * exp( - N * ( ";
            if (con2.node->getNameNum() != datum)
                outFile << "x(" << con2.node->getNameNum() << ")";
            if (con1.node->getNameNum() != datum)
                outFile << " - x(" << con1.node->getNameNum() << ")";
            outFile << " ) )";
            // if (con2.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
            //     JAC[actualName][nameNum2] -= (-1.0 * 1e-16) * (-38.78) / 0.99 * exp(38.78 * (X[con2.node->getNameNum()] - X[con1.node->getNameNum()]));
            // } else if (con2.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
            //     JAC[actualName][nameNum2] -= (-1.0 * 1e-16) * (-38.78) / 0.99 * exp(38.78 * X[con2.node->getNameNum()]);
            // } else {
            //     JAC[actualName][nameNum2] -= (-1.0 * 1e-16) * (-38.78) / 0.99 * exp(38.78 * (-X[con1.node->getNameNum()]));
            // }
            JAC[actualName][nameNum2] -= calculateFe_(X, is, af, n, con2.node->getNameNum(), con1.node->getNameNum(), datum, nameNum2);
        }
        else if ((con1.node->getNameNum() == nameNum1) && (model->getType() == NPN) && (con1.node->getNameNum() == nameNum2)) {
            // double temp1 = (-1.0 * 1e-16) * (-38.78) / 0.99 * exp(38.78 * (X[con2.node->getNameNum()] - X[con1.node->getNameNum()]));
            // double temp2 = (-1.0 * 1e-16) * (-38.78) * exp(38.78 * (X[con0.node->getNameNum()] - X[con1.node->getNameNum()]));
            // double temp3 = (-1.0 * 1e-16) * (-38.78) * exp(38.78 * (X[con2.node->getNameNum()] - X[con1.node->getNameNum()]));
            // double temp4 = (-1.0 * 1e-16) * (-38.78) / 0.5 * exp(38.78 * (X[con0.node->getNameNum()] - X[con1.node->getNameNum()]));
            // JAC[actualName][nameNum2] += temp1 - temp2 - temp3 + temp4;
            outFile << " + IS * N / af * exp( - N * ( ";
            if (con2.node->getNameNum() != datum)
                outFile << "x(" << con2.node->getNameNum() << ")";
            if (con1.node->getNameNum() != datum)
                outFile << " - x(" << con1.node->getNameNum() << ")";
            // if (con2.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
            //     JAC[actualName][nameNum2] += (-1.0 * 1e-16) * (-38.78) / 0.99 * exp(38.78 * (X[con2.node->getNameNum()] - X[con1.node->getNameNum()]));
            // } else if (con2.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
            //     JAC[actualName][nameNum2] += (-1.0 * 1e-16) * (-38.78) / 0.99 * exp(38.78 * X[con2.node->getNameNum()]);
            // } else {
            //     JAC[actualName][nameNum2] += (-1.0 * 1e-16) * (-38.78) / 0.99 * exp(38.78 * (-X[con1.node->getNameNum()]));
            // }
            JAC[actualName][nameNum2] -= calculateFe_(X, is, af, n, con2.node->getNameNum(), con1.node->getNameNum(), datum, nameNum2);
            outFile << " ) )";
            outFile << " - ";
            outFile << "IS * N * exp( - N * ( ";
            if (con0.node->getNameNum() != datum)
                outFile << "x(" << con0.node->getNameNum() << ")";
            if (con1.node->getNameNum() != datum)
                outFile << " - x(" << con1.node->getNameNum() << ")";
            // if (con0.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
            //     JAC[actualName][nameNum2] -= (-1.0 * 1e-16) * (-38.78) * exp(38.78 * (X[con0.node->getNameNum()] - X[con1.node->getNameNum()]));
            // } else if (con0.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
            //     JAC[actualName][nameNum2] -= (-1.0 * 1e-16) * (-38.78) * exp(38.78 * X[con0.node->getNameNum()]);
            // } else {
            //     JAC[actualName][nameNum2] -= (-1.0 * 1e-16) * (-38.78) * exp(38.78 * (-X[con1.node->getNameNum()]));
            // }
            JAC[actualName][nameNum2] += ar * calculateFc_(X, is, ar, n, con0.node->getNameNum(), con1.node->getNameNum(), datum, nameNum2);
            outFile << " ) )";
            outFile << " - ";
            outFile << "IS * N * exp( - N * ( ";
            if (con2.node->getNameNum() != datum)
                outFile << "x(" << con2.node->getNameNum() << ")";
            if (con1.node->getNameNum() != datum)
                outFile << " - x(" << con1.node->getNameNum() << ")";
            // if (con2.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
            //     JAC[actualName][nameNum2] -= (-1.0 * 1e-16) * (-38.78) * exp(38.78 * (X[con2.node->getNameNum()] - X[con1.node->getNameNum()]));
            // } else if (con2.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
            //     JAC[actualName][nameNum2] -= (-1.0 * 1e-16) * (-38.78) * exp(38.78 * X[con2.node->getNameNum()]);
            // } else {
            //     JAC[actualName][nameNum2] -= (-1.0 * 1e-16) * (-38.78) * exp(38.78 * (-X[con1.node->getNameNum()]));
            // }
            JAC[actualName][nameNum2] += af * calculateFe_(X, is, af, n, con2.node->getNameNum(), con1.node->getNameNum(), datum, nameNum2);
            outFile << " ) )";
            outFile << " + ";
            outFile << "IS * N / ar * exp( - N * ( ";
            if (con0.node->getNameNum() != datum)
                outFile << "x(" << con0.node->getNameNum() << ")";
            if (con1.node->getNameNum() != datum)
                outFile << " - x(" << con1.node->getNameNum() << ")";
            outFile << " ) ) ";
            // if (con0.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
            //     JAC[actualName][nameNum2] += (-1.0 * 1e-16) * (-38.78) / 0.5 * exp(38.78 * (X[con0.node->getNameNum()] - X[con1.node->getNameNum()]));
            // } else if (con0.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
            //     JAC[actualName][nameNum2] += (-1.0 * 1e-16) * (-38.78) / 0.5 * exp(38.78 * X[con0.node->getNameNum()]);
            // } else {
            //     JAC[actualName][nameNum2] += (-1.0 * 1e-16) * (-38.78) / 0.5 * exp(38.78 * (-X[con1.node->getNameNum()]));
            // }
            JAC[actualName][nameNum2] -= calculateFc_(X, is, ar, n, con0.node->getNameNum(), con1.node->getNameNum(), datum, nameNum2);
        }
    }
        // TODO：其他情况的处理
        break;

    case VSource:
        // 因为这边是只考虑KCL方程的，所以需要对nameNum2做一个转变【即这里是对未知电流变量求导】
        if (((con0.node->getNameNum() == actualName) && (con1.node->getNameNum() == nameNum2))) {
            JAC[actualName][lastnode + getcompNum()] += 1;
            outFile << " + 1";
        }
        else if((con1.node->getNameNum() == actualName) && (con0.node->getNameNum() == nameNum2)) {
            JAC[actualName][lastnode + getcompNum()] -= 1;
            outFile << " - 1";
        }
        break;

    case ISource:
        break;

    case Diode:
        // TODO：后续学习到再补充处理
        break;

    case Resistor:
        if (((con0.node->getNameNum() == nameNum1) && (con0.node->getNameNum() == nameNum2)) || ((con1.node->getNameNum() == nameNum1) && (con1.node->getNameNum() == nameNum2))) {
            JAC[nameNum1][nameNum2] += 1 / value;
            outFile << " + 1 / " << name;
        } else if (((con0.node->getNameNum() == nameNum1) && (con1.node->getNameNum() == nameNum2)) || ((con1.node->getNameNum() == nameNum1) && (con0.node->getNameNum() == nameNum2))) {
            JAC[nameNum1][nameNum2] -= 1 / value;
            outFile << " - 1 / " << name;
        }
        break;

    case Capacitor:
        // TODO：后续学习到再补充处理
        break;

    case Inductor:
        cerr << "This section is not completed" << endl
             << "PROGRAM ENDED ABNORMALLY!" << endl;
        exit(0);
        break;
    }
    return;
}

/****************************************************************************************************************************************************************/

/****************************************************************************************************************************************************************/
/****************************************************************************************************************************************************************/

Node::Node(int Num) {
    next = NULL;
    nodeNum = Num;
    conCount = 0;
    conList = NULL;
    nameNum = NA;
}

Node::~Node() {};

// 获取节点在链上的编号【其实用于判断节点号一样的时候，用这个namenum也行】
int Node::getNum() {
    return nodeNum;
}

// 设置节点的编号
void Node::setNameNum(int numIn) {
    nameNum = numIn;
}

// 获取节点在网表中命名的编号
int Node::getNameNum() {
    return nameNum;
}

// 获取本节点连接的器件数目
int Node::getCount() {
    return conCount;
}

// 获取本节点相关的连接关系
Connections *Node::getConList() {
    return conList;
}

// 将组件连接到节点上
void Node::connect(int conNumIn, Component *compIn) {
    Connections *conPtr;
    conCount++;
    if (conList == NULL) { //如果当前节点的连接组件链表为空
        conList = new Connections;
        conList->next = NULL;
        conList->conNum = conNumIn; //这个连接关系的端口是conNumIn号端口
        conList->comp = compIn;     //指向了传过来的组件
    }
    else {
        conPtr = conList;
        while (conPtr->next != NULL)
            conPtr = conPtr->next;
        conPtr->next = new Connections;
        conPtr = conPtr->next;
        conPtr->next = NULL;
        conPtr->conNum = conNumIn;
        conPtr->comp = compIn;
    }
}

// 获取下一个节点指针
Node *Node::getNext() {
    return next;
}

// 添加节点到链上，并完善连接关系
void Node::setNext(Node *nodeIn) {
    next = nodeIn;
}

// 打印节点信息
void Node::printMessage(ofstream &outFile){
    Connections *conList = getConList();
    while(conList != NULL){
        conList->comp->printMessage(outFile, conList->conNum);
        conList = conList->next;
    }
}

/***********************************************************************************************************************************/

/*
 * 生成本节点的方程（KCL）
 * @param outFile 输出文件流
 * @param F_x F(X)方程组
 * @param X 待求解向量
 * @param datum 接地节点编号
 * @param lastnode 最大节点编号
 * @param isMNA 是否是MNA节点
 */
void Node::genKCLEquation(ofstream &outFile, vector<double> &F_x, vector<double>& X, int datum, int lastnode, Boolean isMNA) {
    int VSourceID = NA;
    Connections *conList = getConList();
    // TODO：这里假设了本节点只连接了一个电压源来处理的，后续若有问题再做处理
    while(conList != NULL) {
        if(conList->comp->getType() == VSource){
            VSourceID = conList->comp->getcompNum();
            break;
        }
        conList = conList->next;
    }
    if(isMNA == TRUE) {
        outFile << "F(" << lastnode + VSourceID << ") = ";
    } else{
        outFile << "F(" << nameNum << ") = ";
    }
    conList = getConList();
    while(conList != NULL){
        conList->comp->genKCLEquation(outFile, F_x, X, datum, lastnode, nameNum, isMNA == TRUE ? lastnode + VSourceID : NA);
        conList = conList->next;
    }
    outFile << endl;
}

/*
 * 生成本组件相关的方程（KCL）的导数
 * @param outFile 输出文件流
 * @param JAC F(X)方程组的雅可比矩阵
 * @param X 待求解向量
 * @param nameNum2 二维编号
 * @param datum 接地节点编号
 * @param lastnode 最大节点编号
 * @param isMNA 是否是MNA节点（一维）
 */
void Node::genKCLJAC(ofstream &outFile, vector<vector<double>> &JAC, vector<double>& X, int nameNum2, int datum, int lastnode, Boolean isMNA) {
    int VSourceID = NA;
    Boolean isAnotherCon = FALSE;
    Connections *conList = getConList();
    // TODO：这里假设了本节点只连接了一个电压源来处理的，后续若有问题再做处理
    while(conList != NULL) {
        if(conList->comp->getType() == VSource){
            VSourceID = conList->comp->getcompNum();
            isAnotherCon = (conList->comp->getConVal(0) == nameNum2 || conList->comp->getConVal(1) == nameNum2) ? TRUE : FALSE;
            break;
        }
        conList = conList->next;
    }
    if(isMNA == TRUE) {
        outFile << "JAC(" << lastnode + VSourceID << "," << nameNum2 << ") = ";
    } else { // 如果是两个均非接地的节点之间的电压源，那么其中一个作为KVL方程后，会有一个位置电流变量，处理如下：
        if(isAnotherCon == TRUE && nameNum2 != nameNum)
            outFile << "JAC(" << nameNum << "," << VSourceID + lastnode << ") = ";
        else
            outFile << "JAC(" << nameNum << "," << nameNum2 << ") = ";
    }
    conList = getConList();
    while(conList != NULL){
        conList->comp->genKCLJAC(outFile, JAC, X, nameNum, nameNum2, datum, lastnode, isMNA == TRUE ? lastnode + VSourceID : NA);
        conList = conList->next;
    }
    outFile << endl;
}

NodeHead::NodeHead() {
    nodeList = NULL;
    nodeCount = 0;
}

NodeHead::~NodeHead(){};

// 添加节点到链上
Node *NodeHead::addNode() {
    Node *nodePtr;
    nodeCount++;
    if (nodeList == NULL) {
        nodeList = new Node(nodeCount);
        return nodeList;
    }
    else {
        nodePtr = nodeList;
        while (nodePtr->getNext() != NULL)
            nodePtr = nodePtr->getNext();
        nodePtr->setNext(new Node(nodeCount));
        return nodePtr->getNext();
    }
}

// 获取节点链表上的节点总数
int NodeHead::getCount() {
    return nodeCount;
}

// 根据节点编号获取节点实体
Node *NodeHead::getNode(int nodeNum) {
    Node *nodePtr = nodeList;
    // need check that nodeNum does not exceed node count
    // 因为在链上节点是按照顺序来的，逐个依序编号的
    for (int a = 0; a < nodeNum; a++)
        nodePtr = nodePtr->getNext();
    return nodePtr;
}

CompHead::CompHead() {
    compList = NULL;
    mCount = 0;
    bCount = 0;
    iCount = 0;
    rCount = 0;
    dCount = 0;
    cCount = 0;
    vSCount = 0;
    iSCount = 0;
}

CompHead::~CompHead(){};

// 添加器件到链表
void CompHead::addComp(Component *component) {
    Component *compPtr;
    switch (component->getType()) {
    case ISource:
        iSCount++;
        if (compList == NULL) {
            compList = component;
            component->setNum(iSCount);
        }
        else {
            compPtr = compList;
            while (compPtr->getNext() != NULL)
                compPtr = compPtr->getNext();
            compPtr->setNext(component);
            compPtr = compPtr->getNext();
            compPtr->setNum(iSCount);
        }
        break;
    case VSource:
        vSCount++;
        if (compList == NULL) {
            compList = component;
            component->setNum(vSCount);
        }
        else {
            compPtr = compList;
            while (compPtr->getNext() != NULL)
                compPtr = compPtr->getNext();
            compPtr->setNext(component);
            compPtr = compPtr->getNext();
            compPtr->setNum(vSCount);
        }
        break;
    case Resistor:
        rCount++;
        if (compList == NULL)
        {
            compList = component;
            component->setNum(rCount);
        }
        else {
            compPtr = compList;
            while (compPtr->getNext() != NULL)
                compPtr = compPtr->getNext();
            compPtr->setNext(component);
            compPtr = compPtr->getNext();
            compPtr->setNum(rCount);
        }
        break;
    case MOSFET:
        mCount++;
        if (compList == NULL) {
            compList = component;
            component->setNum(mCount);
        }
        else {
            compPtr = compList;
            while (compPtr->getNext() != NULL)
                compPtr = compPtr->getNext();
            compPtr->setNext(component);
            compPtr = compPtr->getNext();
            compPtr->setNum(mCount);
        }
        break;
    case BJT:
        bCount++;
        if (compList == NULL) {
            compList = component;
            component->setNum(bCount);
        }
        else {
            compPtr = compList;
            while (compPtr->getNext() != NULL)
                compPtr = compPtr->getNext();
            compPtr->setNext(component);
            compPtr = compPtr->getNext();
            compPtr->setNum(bCount);
        }
        break;
    case Diode:
        dCount++;
        if (compList == NULL) {
            compList = component;
            component->setNum(dCount);
        }
        else {
            compPtr = compList;
            while (compPtr->getNext() != NULL)
                compPtr = compPtr->getNext();
            compPtr->setNext(component);
            compPtr = compPtr->getNext();
            compPtr->setNum(dCount);
        }
        break;
    case Capacitor:
        cCount++;
        if (compList == NULL) {
            compList = component;
            component->setNum(cCount);
        }
        else {
            compPtr = compList;
            while (compPtr->getNext() != NULL)
                compPtr = compPtr->getNext();
            compPtr->setNext(component);
            compPtr = compPtr->getNext();
            compPtr->setNum(cCount);
        }
        break;
    case Inductor:
        iCount++;
        if (compList == NULL) {
            compList = component;
            component->setNum(iCount);
        }
        else {
            compPtr = compList;
            while (compPtr->getNext() != NULL)
                compPtr = compPtr->getNext();
            compPtr->setNext(component);
            compPtr = compPtr->getNext();
            compPtr->setNum(iCount);
        }
        break;
    }
}

// 获取传入类型的器件的总数
int CompHead::getCount(CompType type) {
    switch (type) {
    case ISource:
        return iSCount;
    case VSource:
        return vSCount;
    case Resistor:
        return rCount;
    case Diode:
        return dCount;
    case MOSFET:
        return mCount;
    case BJT:
        return bCount;
    case Capacitor:
        return cCount;
    case Inductor:
        return iCount;
    }
    return NA;
}

// 依据器件编号（链上位置），获取器件实体
Component *CompHead::getComp(int compNum) {
    Component *compPtr = compList;
    for (int a = 0; a < compNum; a++)
        compPtr = compPtr->getNext();
    return compPtr;
}

Model::Model(char *nameIn, TranType typeIn, double isIn, double bfIn, double brIn, double tempIn) {
    // 由于网表参数的正负值问题，可能会出现结果的错误。
    // 本次设计是ishen为负值
    strcpy(name, nameIn);
    type = typeIn;
    is = isIn;
    br = brIn;
    bf = bfIn;
    af = bfIn / (bfIn + 1);
    ar = brIn / (brIn + 1);
    temp = tempIn;
    next = NULL;
}

Model::~Model() {}

// 获取model类型
TranType Model::getType() {
    return type;
}

// 获取模型名称
char *Model::getName() {
    return name;
}

// 获取模型is参数值
double Model::getIs() {
    return is;
}

// 获取模型bf参数值
double Model::getBf() {
    return bf;
}

// 获取模型br参数值
double Model::getBr() {
    return br;
}

// 获取模型af参数值
double Model::getAf() {
    return af;
}

// 获取模型ar参数值
double Model::getAr() {
    return ar;
}

// 获取模型temp参数值
double Model::getTemp() {
    return temp;
}

// 获取N值【注意这里设定了一个恒定的N值，如果temp为默认值NA的话，那么采用38.78】
double Model::getN() {
    return temp == NA ? 38.78 : (Q / (K * getTemp()));
}

// 添加下一项
void Model::setNext(Model *nextIn) {
    next = nextIn;
}

// 获取下一项指针
Model *Model::getNext() {
    return next;
}

// 模型链表初始化
ModelHead::ModelHead() {
    modelList = NULL;
}

// 添加模型到链上
void ModelHead::addModel(Model *modelIn) {
    Model *modelPtr;
    if (modelList == NULL) {
        modelList = modelIn;
    }
    else {
        modelPtr = modelList;
        while (modelPtr->getNext() != NULL) {
            modelPtr = modelPtr->getNext();
        }
        modelPtr->setNext(modelIn);
    }
}

// 依据传入名称获取模型实体
Model *ModelHead::getModel(char *nameIn) {
    Model *modelPtr = modelList;
    while (strcmp(modelPtr->getName(), nameIn)) {
        modelPtr = modelPtr->getNext();
    }
    return modelPtr;
}