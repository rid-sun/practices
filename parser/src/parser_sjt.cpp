//  PARSER.CC  Ver. 1.0
//  Program to extract Nodal equations from Spice Netlist.  Ed Chan

#include "parser_sjt.h"

void generateMatrix(NodeHead &nodeList, CompHead &compList, ModelHead &modelList, vector<double> &F_x, vector<double> &X, vector<vector<double>> &JAC, string &outFileName, int datum, int lastnode, int step) {
    string outFile_name;
    ofstream outFile;
    Node *nodePtr, *nodePtr1, *nodePtr2;
    Component *comp;

    // 打开输出文件流
    outFile_name = (outFileName + "_equation") + to_string(step);
    outFile.open(outFile_name + ".txt", ios::out);
    if(!outFile){
        cerr << "打开文件出现错误，异常退出！" << endl;
        exit(1);
    }

    // 1. 按节点分析法生成节点方程，如有电压源，则引入补偿方程等，即转为MNA方程
    nodePtr = nodeList.getNode(0);
    while (nodePtr != NULL) {
        // Connections *conList = nodePtr->getConList();
        // Boolean isHaveVSource = FALSE;
        // while(conList != NULL) {
        //     // 按照Eric的代码，这里也是有电感器存在的判断
        //     // TODO：确认是否需要加电感器的判断
        //     if(conList->comp->getType() == VSource) {
        //         isHaveVSource = TRUE;
        //         break;
        //     }
        //     conList = conList->next;
        // }
        // // 对于电压源连的节点，另外处理
        // if ( isHaveVSource == FALSE) { 
        //     nodePtr->genNodalEquation(outFile, F_x, X, datum, lastnode);
        // } else {
        //     nodePtr->genSpecialNodal(outFile, F_x, X, datum, lastnode);
        // }
        if(nodePtr->getNameNum() != datum) {
            nodePtr->genNodalEquation(outFile, F_x, X, datum, lastnode);
        }
        nodePtr = nodePtr->getNext();
    }
    outFile.close();
    outFile.clear();

    
    
    // // 2. 仍然扫描节点链表，找出连接电压源的额节点
    // comp = compList.getComp(0);
    // while (comp != NULL) {
    //     if(comp->getType() != VSource) continue;
    //     comp->genMNAEquation(outFile, F_x, X, datum, lastnode);
    // }
    // nodePtr = nodeList.getNode(0);
    // 3. 对按节点分析法生成的节点方程求导
    // 4. 扫描源，若有则对引入改进节点方程部分求导



    
    // 打开输出文件流
    outFile_name = (outFileName + "_jacobian") + to_string(step);
    outFile.open(outFile_name + ".txt", ios::out);
    if(!outFile){
        cerr << "打开文件出现错误，异常退出！" << endl;
        exit(1);
    }

    // 2. 生成上面求得方程的雅可比矩阵
    nodePtr1 = nodeList.getNode(0);
    while(nodePtr1 != NULL) {
        if(nodePtr1->getNameNum() != datum) {
            nodePtr2 = nodeList.getNode(0);
            while(nodePtr2 != NULL) {
                if(nodePtr2->getNameNum() != datum)
                    nodePtr1->genNodalJAC(outFile, JAC, X, nodePtr2, datum, lastnode);
                nodePtr2 = nodePtr2->getNext();
            }
        }
        nodePtr1 = nodePtr1->getNext();
    }
    outFile.close();
    outFile.clear();

    cout << "==================The solution of the "<<step<<"th iteration is completed.==================" << endl;
}

void parseNetList(NodeHead &nodeList, CompHead &compList, ModelHead &modelList, int& datum, int& lastnode, string& inFileName, string& outFileName) {
    ifstream inFile;
    ofstream outFile;

    // 分析过程中所需的变量
    char buf[BufLength], buf1[BufLength], buf2[BufLength], buf3[BufLength], nameBuf[NameLength];
    char *bufPtr, *charPtr1, *charPtr2;
    int intBuf1, intBuf2, intBuf3, intBuf4, eqNum = NA;
    double douBuf1, douBuf2, douBuf3, douBuf4;
    CompType typeBuf;
    Component *compPtr, *compPtr1, *compPtr2;
    Node *nodePtr, *nodePtr1, *nodePtr2;
    Model *modelPtr;
    TranType TtypeBuf;
    // EquaType eqType = Modified;

    // 1. 方程类型选择 【这部分暂且忽略，默认MNA的方程】
    // if (eqNum == NA) {
    //     while ((eqNum != 1) && (eqNum != 2)) {
    //         cout << "Available Equations Types Are:" << endl
    //              << " <1>  Nodal" << endl
    //              << " <2>  Modified Nodal" << endl
    //              << "Please enter your choice <1, 2>:" << endl;
    //         cin >> buf;
    //         eqNum = atoi(buf);
    //     }
    //     if (eqNum == 1)
    //         eqType = Nodal;
    //     else if (eqNum == 2)
    //         eqType = Modified;
    // }

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
            while (charPtr1 != NULL) {    //若有参数
                // 下面处理就是四种参数的值，分别对应就行。stripString是将字符串转成浮点数double的
                if ((charPtr1[0] == 'I') && (charPtr1[1] == 'S') && (charPtr1[2] == '=')) {
                    douBuf1 = stripString(charPtr1);
                }
                if ((charPtr1[0] == 'B') && (charPtr1[1] == 'F') && (charPtr1[2] == '=')) {
                    douBuf2 = stripString(charPtr1);
                }
                if ((charPtr1[0] == 'B') && (charPtr1[1] == 'R') && (charPtr1[2] == '=')) {
                    douBuf3 = stripString(charPtr1);
                }
                if ((charPtr1[0] == 'T') && (charPtr1[1] == 'E') && (charPtr1[4] == '=')) {
                    douBuf4 = stripString(charPtr1);
                }
                charPtr1 = strtok(NULL, " ");
            }
            modelPtr = new Model(buf2, TtypeBuf, douBuf1, douBuf2, douBuf3, douBuf4);
            modelList.addModel(modelPtr);
        }
        inFile.getline(buf, BufLength);
    } //模型已经扫描处理完毕
    inFile.close();
    inFile.clear();// 这点注意一下
    inFile.open(inFileName, ios::in); //现在开始扫描处理其他的信息

    // 4.2 器件遍历扫描，并入链
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

    /************************************************************/
    /** 对于起始节点编号的确定是一个问题。常见的是接地点为0号节点， **/
    /** 也即为第一个节点。但似乎有约定是第一个源的第二个端口所连节  **/
    /** 点作为起始节点？有待考证。或者个人指定起始节点             **/
    /************************************************************/

    // 4.3 依据器件链表，进行节点遍历扫描，节点成链并与器件建立联系
    compPtr1 = compList.getComp(0);
    while (compPtr1 != NULL) {
        for (int b = 0; b < 3; b++) {
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

    cout << "=====================The netlist is parsed completely=========================" << endl;

    return;
}

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
};

// Print the linked list of components to check
void printComponents(Component *compPtr) {
    char compTypeName[6];
    cout << endl
         << "Components: " << endl
         << endl;
    while (compPtr != NULL) {
        strcpy(compTypeName, strComponentType(compPtr));
        cout << "->" << compTypeName << compPtr->getcompNum();
        compPtr = compPtr->getNext();
    }
    cout << endl;
    return;
}

void printNodes(Node *nodePtr, int compFlag) {
    Connections *conPtr;
    cout << endl
         << "Nodes: " << endl
         << endl;
    while (nodePtr != NULL) {
        if (compFlag == 0) { // It is printed just the names of the nodes
            cout << "-> " << nodePtr->getNameNum();
        } else if (compFlag == 1) { // It is printed the nodes and the connections
            cout << "-> " << nodePtr->getNameNum() << " {";
            conPtr = nodePtr->getConList();
            while (conPtr->next != NULL) {
                cout << strComponentType(conPtr->comp) << conPtr->comp->getcompNum() << ", ";
                conPtr = conPtr->next;
            }
            cout << strComponentType(conPtr->comp) << conPtr->comp->getcompNum() << '}' << endl;
        } else {
            cout << "Invalid value for compFlag. (0) to print just nodes, (1) to print nodes and connections!";
            exit(1);
        }
        nodePtr = nodePtr->getNext();
    }
    return;
}

char *strComponentType(Component *compPtr) {
    char *compTypeName = new char[6];
    switch (compPtr->getType()) {
    case VSource:
        strcpy(compTypeName, "V");
        break;
    case Resistor:
        strcpy(compTypeName, "R");
        break;
    case BJT:
        strcpy(compTypeName, "T");
        break;
    case MOSFET:
        strcpy(compTypeName, "M");
        break;
    case ISource:
        strcpy(compTypeName, "I");
        break;
    case Inductor:
        strcpy(compTypeName, "ind");
        break;
    case Diode:
        strcpy(compTypeName, "Diode");
        break;
    case Capacitor:
        strcpy(compTypeName, "Cap");
        break;
    }
    return compTypeName;
}

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

CompType Component::getType() {
    return type;
}

int Component::getNum() {
    return compNum;
}

int Component::getcompNum() {
    return compNum;
}

Component *Component::getNext() {
    return next;
}

double Component::getVal() {
    return value;
}

void Component::setNext(Component *nextIn) {
    next = nextIn;
}

void Component::setNum(int numIn) {
    compNum = numIn;
}

/*
 * 获得该器件某端口连接的节点号
 * @param conNum “端口”号
 */
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

/*
 * 判断端口是否可用
 * @param conNum “端口”号
 */
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
        outFile << "        value: IS = " << model->getIs() << "  BF = " << model->getBf() << "  BR = " << model->getBr() << "  N = " << (Q / (K * model->getTemp())) << endl;
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

void Component::genEquation(ofstream &outFile, vector<double> &F_x, vector<double>& X, int datum, int lastnode, int nameNum) {
    switch (type) {
    case MOSFET:
        //这里存疑？是否MOSFET与BJT是一个特性的东西
    case BJT:
        if ((con0.node->getNameNum() == nameNum) && (model->getType() == NPN)) {
            double fc = (-1.0 * 1e-16) / 0.5 * (exp(-38.78 * (X[con0.node->getNameNum()] - X[con1.node->getNameNum()])) - 1);
            double fe = (-1.0 * 1e-16) / 0.5 * (exp(-38.78 * (X[con2.node->getNameNum()] - X[con1.node->getNameNum()])) - 1);
            F_x[nameNum] = F_x[nameNum] + fc - 0.99 * fe;
            outFile << " + " << name << "_Ic";
        }
        if ((con2.node->getNameNum() == nameNum) && (model->getType() == NPN)) {
            double fc = (-1.0 * 1e-16) / 0.5 * (exp(-38.78 * (X[con0.node->getNameNum()] - X[con1.node->getNameNum()])) - 1);
            double fe = (-1.0 * 1e-16) / 0.5 * (exp(-38.78 * (X[con2.node->getNameNum()] - X[con1.node->getNameNum()])) - 1);
            F_x[nameNum] = F_x[nameNum] + fe - 0.5 * fc;
            outFile << " + " << name << "_Ie";     
        }
        // if((con1.node->getNameNum() == nameNum) && (model->getType() == NPN)) {
        //     Connections *conList = con1.node->getConList();
        //     Boolean isOnlyBJT = TRUE;
        //     while(conList != NULL){
        //         // 这边需要进行一个判定，判定节点还有连接电阻之类的器件
        //         // Todo: 是不是只有这种情况呢？需要后续了解处理
        //         if(conList->comp->getType() == Resistor) {
        //             isOnlyBJT = FALSE;
        //             break;
        //         }
        //         conList = conList->next;
        //     }
        //     if(isOnlyBJT == FALSE){
        //         double fc = (-1.0 * 1e-16) / 0.5 * (exp(-38.78 * (X[con0.node->getNameNum()] - X[con1.node->getNameNum()])) - 1);
        //         double fe = (-1.0 * 1e-16) / 0.5 * (exp(-38.78 * (X[con2.node->getNameNum()] - X[con1.node->getNameNum()])) - 1);
        //         F_x[nameNum] -= fc - 0.99 * fe + fe - 0.5 * fc;
        //         outFile << " - " << name << "_Ie"
        //                 << " - " << name << "_Ic";
        //     }
        // }
        if((con1.node->getNameNum() == nameNum) && (model->getType() == NPN)) {
            double fc = (-1.0 * 1e-16) / 0.5 * (exp(-38.78 * (X[con0.node->getNameNum()] - X[con1.node->getNameNum()])) - 1);
            double fe = (-1.0 * 1e-16) / 0.5 * (exp(-38.78 * (X[con2.node->getNameNum()] - X[con1.node->getNameNum()])) - 1);
            F_x[nameNum] -= fc - 0.99 * fe + fe - 0.5 * fc;
            outFile << " - " << name << "_Ie"
                    << " - " << name << "_Ic";
        }
        // Todo: 其他情况是知识完善后，再做补充处理
        break;

    case VSource://KVL方程
        if (con0.node->getNameNum() != datum && con1.node->getNameNum() != datum){
            // TODO：似乎这种情况属于eric所说的float source？
        } else {
            // 因为datum节点是不会生成方程从而进入这里的，所以可以直接像下面这样写
            if (con0.node->getNameNum() == nameNum) {
                F_x[nameNum] += X[nameNum] - value;
                outFile << " + "
                        << "x(" << nameNum << ") - " << name;
            } else if (con1.node->getNameNum() == nameNum) {
                F_x[nameNum] += X[nameNum] + value;
                outFile << " + "
                        << "x(" << nameNum << ") + " << name;
            }
        }
        break;

    case ISource:
        if (con0.node->getNameNum() == nameNum) {
            F_x[nameNum] += value;
            outFile << " + " << name;
        } else if (con1.node->getNameNum() == nameNum) {
            F_x[nameNum] -= value;
            outFile << "-" << name;
        }
        break;

    case Diode:
        // if (con0.node->getNameNum() == nameNum)
        // {
        //     outFile << " (" << name << "IS "
        //             << ")*(exp("
        //             << name << "N*(";
        //     if (con0.node->getNameNum() != datum)
        //         outFile << "X(" << con0.node->getNameNum() << ')';
        //     if (con1.node->getNameNum() != datum)
        //         outFile << "-X(" << con1.node->getNameNum() << ')';
        //     outFile << ")) -1) ";
        // }
        // else if (con1.node->getNameNum() == nameNum)
        // {
        //     outFile << " (-" << name << "IS "
        //             << ")*(exp("
        //             << name << "N*(";
        //     if (con0.node->getNameNum() != datum)
        //         outFile << "X(" << con0.node->getNameNum() << ')';
        //     if (con1.node->getNameNum() != datum)
        //         outFile << "-X(" << con1.node->getNameNum() << ')';
        //     outFile << ")) -1) ";
        // }
        break;

    case Resistor:
        if (con0.node->getNameNum() == nameNum) {
            F_x[nameNum] += (X[nameNum] - X[con1.node->getNameNum()]) / value;
            outFile << " + (";
            if (con0.node->getNameNum() != datum)
                outFile << " x(" << con0.node->getNameNum() << ')';
            if (con1.node->getNameNum() != datum)
                outFile << " - x(" << con1.node->getNameNum() << ')';
            outFile << " ) / " << name;
        }
        if (con1.node->getNameNum() == nameNum) {
            F_x[nameNum] += (X[nameNum] - X[con0.node->getNameNum()]) / value;
            outFile << " + (";
            if (con1.node->getNameNum() != datum)
                outFile << " x(" << con1.node->getNameNum() << ')';
            if (con0.node->getNameNum() != datum)
                outFile << " - x(" << con0.node->getNameNum() << ')';
            outFile << " ) / " << name;
        }
        break;

    case Capacitor:
        // outFile << " 0 ";
        break;

    case Inductor:
        // if (con0.node->getNum() == nodeNum)
        //     outFile << " Il" << compNum << " ";
        // else if (con1.node->getNum() == nodeNum)
        break;
    }
    return;
}

// void Component::genMNAEquation(ofstream &outFile, vector<double> &F_x, vector<double> X, int datum, int lastnode) {
//     // assert：假定电压源的一端是datum，也即接地端
//     // TODO：可能会有两端不接地的情况？也就是eric代码中的supernode？？？
//     outFile << "F(" << getcompNum() + lastnode << ") = x(" << getcompNum() + lastnode << ")";
//     if(con0.node->getNameNum() != datum) {
//         genEquation(outFile, F_x, X, datum, lastnode, con0.node->getNameNum());
//         outFile << endl;
//     } else if(con0.node->getNameNum() != datum) {
//         genEquation(outFile, F_x, X, datum, lastnode, con0.node->getNameNum());
//         outFile << endl;   
//     }
// }

void Component::genJAC(ofstream &outFile, vector<vector<double>> &JAC, vector<double>& X, int nameNum1, int nameNum2, int datum, int lastnode) {
    switch (type) {
    case MOSFET:
        std::cout << "genJAC for MOSFETs not implemented" << std::endl
                  << "PROGRAM ENDED ABNORMALLY!" << std::endl;
        exit(0);
        break;

    case BJT:
        if ((con0.node->getNameNum() == nameNum1) && (model->getType() == NPN) && (con0.node->getNameNum() == nameNum2)) {
            JAC[nameNum1][nameNum2] += (-1.0 * 1e-16) * (-38.78) / 0.5 * exp(38.78 * (X[con0.node->getNameNum()] - X[con1.node->getNameNum()]));
            outFile << " + IS * N / ar * exp( - N * ( ";
            if (con0.node->getNameNum() != datum)
                outFile << "x(" << con0.node->getNameNum() << ")";
            if (con1.node->getNameNum() != datum)
                outFile << " - x(" << con1.node->getNameNum() << ")";
            outFile << " ) )";
        }
        else if ((con0.node->getNameNum() == nameNum1) && (model->getType() == NPN) && (con2.node->getNameNum() == nameNum2)) {
            JAC[nameNum1][nameNum2] += (1.0 * 1e-16) * (-38.78) * exp(38.78 * (X[con2.node->getNameNum()] - X[con1.node->getNameNum()]));
            outFile << " - IS * N * exp( - N * ( ";
            if (con2.node->getNameNum() != datum)
                outFile << "x(" << con2.node->getNameNum() << ")";
            if (con1.node->getNameNum() != datum)
                outFile << " - x(" << con1.node->getNameNum() << ")";
            outFile << " ) )";
        }
        else if ((con0.node->getNameNum() == nameNum1) && (model->getType() == NPN) && (con1.node->getNameNum() == nameNum2)) {
            double temp1 = (-1.0 * 1e-16) * (-38.78) * exp(38.78 * (X[con2.node->getNameNum()] - X[con1.node->getNameNum()]));
            double temp2 = (-1.0 * 1e-16) * (-38.78) / 0.5 * exp(38.78 * (X[con0.node->getNameNum()] - X[con1.node->getNameNum()]));
            JAC[nameNum1][nameNum2] += temp1 - temp2;
            outFile << " + IS * N * exp( - N * ( ";
            if (con2.node->getNameNum() != datum)
                outFile << "x(" << con2.node->getNameNum() << ")";
            if (con1.node->getNameNum() != datum)
                outFile << " - x(" << con1.node->getNameNum() << ")";
            outFile << " ) )";
            outFile << " - ";
            outFile << "IS * N / ar * exp( - N * ( ";
            if (con0.node->getNameNum() != datum)
                outFile << "x(" << con0.node->getNameNum() << ")";
            if (con1.node->getNameNum() != datum)
                outFile << " - x(" << con1.node->getNameNum() << ")";
            outFile << " ) )";
        }
        else if ((con2.node->getNameNum() == nameNum1) && (model->getType() == NPN) && (con0.node->getNameNum() == nameNum2)) {
            JAC[nameNum1][nameNum2] -= (-1.0 * 1e-16) * (-38.78) * exp(38.78 * (X[con0.node->getNameNum()] - X[con1.node->getNameNum()]));
            outFile << " - IS * N * exp( - N * ( ";
            if (con0.node->getNameNum() != datum)
                outFile << "x(" << con0.node->getNameNum() << ')';
            if (con1.node->getNameNum() != datum)
                outFile << " - x(" << con1.node->getNameNum() << ')';
            outFile << " ) )";
        }
        else if ((con2.node->getNameNum() == nameNum1) && (model->getType() == NPN) && (con2.node->getNameNum() == nameNum2)) {
            JAC[nameNum1][nameNum2] += (-1.0 * 1e-16) * (-38.78) / 0.99 * exp(38.78 * (X[con2.node->getNameNum()] - X[con1.node->getNameNum()]));
            outFile << "IS * N / af * exp( - N * ( ";
            if (con2.node->getNameNum() != datum)
                outFile << "x(" << con2.node->getNameNum() << ')';
            if (con1.node->getNameNum() != datum)
                outFile << " - x(" << con1.node->getNameNum() << ')';
            outFile << " ) )";
        }
        else if ((con2.node->getNameNum() == nameNum1) && (model->getType() == NPN) && (con1.node->getNameNum() == nameNum2)) {
            double temp2 = (-1.0 * 1e-16) * (-38.78) * exp(38.78 * (X[con0.node->getNameNum()] - X[con1.node->getNameNum()]));
            double temp1 = (-1.0 * 1e-16) * (-38.78) / 0.99 * exp(38.78 * (X[con2.node->getNameNum()] - X[con1.node->getNameNum()]));
            JAC[nameNum1][nameNum2] += temp2 - temp1;
            outFile << " - IS * N / af * exp( - N * ( ";
            if (con2.node->getNameNum() != datum)
                outFile << "x(" << con2.node->getNameNum() << ')';
            if (con1.node->getNameNum() != datum)
                outFile << "-x(" << con1.node->getNameNum() << ')';
            outFile << " ) )";
            outFile << " + ";
            outFile << "IS * N * exp( - N * ( ";
            if (con0.node->getNameNum() != datum)
                outFile << "x(" << con0.node->getNameNum() << ')';
            if (con1.node->getNameNum() != datum)
                outFile << " - x(" << con1.node->getNameNum() << ')';
            outFile << " ) )";
        }
        else if ((con1.node->getNameNum() == nameNum1) && (model->getType() == NPN) && (con0.node->getNameNum() == nameNum2)) {
            double temp1 = (-1.0 * 1e-16) * (-38.78) * exp(38.78 * (X[con0.node->getNameNum()] - X[con1.node->getNameNum()]));
            double temp2 = (-1.0 * 1e-16) * (-38.78) / 0.5 * exp(38.78 * (X[con0.node->getNameNum()] - X[con1.node->getNameNum()]));
            JAC[nameNum1][nameNum2] += temp1 - temp2;
            outFile << " + IS * N * exp( - N * ( ";
            if (con0.node->getNameNum() != datum)
                outFile << "x(" << con0.node->getNameNum() << ')';
            if (con1.node->getNameNum() != datum)
                outFile << " - x(" << con1.node->getNameNum() << ')';
            outFile << " + IS * N / ar * exp( - N * ( ";
            if (con0.node->getNameNum() != datum)
                outFile << "x(" << con0.node->getNameNum() << ')';
            if (con1.node->getNameNum() != datum)
                outFile << " - x(" << con1.node->getNameNum() << ')';
            outFile << " ) ) ";
        }
        else if ((con1.node->getNameNum() == nameNum1) && (model->getType() == NPN) && (con2.node->getNameNum() == nameNum2)) {
            double temp1 = (-1.0 * 1e-16) * (-38.78) / 0.99 * exp(38.78 * (X[con2.node->getNameNum()] - X[con1.node->getNameNum()]));
            double temp2 = (-1.0 * 1e-16) * (-38.78) * exp(38.78 * (X[con2.node->getNameNum()] - X[con1.node->getNameNum()]));
            JAC[nameNum1][nameNum2] += temp2 - temp1;
            outFile << " + IS * N * exp( - N * ( ";
            if (con2.node->getNameNum() != datum)
                outFile << "x(" << con2.node->getNameNum() << ')';
            if (con1.node->getNameNum() != datum)
                outFile << " - x(" << con1.node->getNameNum() << ')';
            outFile << " ) )";
            outFile << " - ";
            outFile << "IS * N / af * exp( - N * ( ";
            if (con2.node->getNameNum() != datum)
                outFile << "x(" << con2.node->getNameNum() << ')';
            if (con1.node->getNameNum() != datum)
                outFile << " - x(" << con1.node->getNameNum() << ')';
            outFile << " ) )";
        }
        else if ((con1.node->getNameNum() == nameNum1) && (model->getType() == NPN) && (con1.node->getNameNum() == nameNum2)) {
            double temp1 = (-1.0 * 1e-16) * (-38.78) / 0.99 * exp(38.78 * (X[con2.node->getNameNum()] - X[con1.node->getNameNum()]));
            double temp2 = (-1.0 * 1e-16) * (-38.78) * exp(38.78 * (X[con0.node->getNameNum()] - X[con1.node->getNameNum()]));
            double temp3 = (-1.0 * 1e-16) * (-38.78) * exp(38.78 * (X[con2.node->getNameNum()] - X[con1.node->getNameNum()]));
            double temp4 = (-1.0 * 1e-16) * (-38.78) / 0.5 * exp(38.78 * (X[con0.node->getNameNum()] - X[con1.node->getNameNum()]));
            JAC[nameNum1][nameNum2] += temp1 - temp2 - temp3 + temp4;
            outFile << " + IS * N / af * exp( - N * ( ";
            if (con2.node->getNameNum() != datum)
                outFile << "x(" << con2.node->getNameNum() << ')';
            if (con1.node->getNameNum() != datum)
                outFile << " - x(" << con1.node->getNameNum() << ')';
            outFile << " ) )";
            outFile << " - ";
            outFile << "IS * N * exp( - N * ( ";
            if (con0.node->getNameNum() != datum)
                outFile << "x(" << con0.node->getNameNum() << ')';
            if (con1.node->getNameNum() != datum)
                outFile << "-x(" << con1.node->getNameNum() << ')';
            outFile << " ) )";
            outFile << " - ";
            outFile << "IS * N * exp( - N * ( ";
            if (con2.node->getNameNum() != datum)
                outFile << "x(" << con2.node->getNameNum() << ')';
            if (con1.node->getNameNum() != datum)
                outFile << " - x(" << con1.node->getNameNum() << ')';
            outFile << " ) )";
            outFile << " + ";
            outFile << "IS * N / ar * exp( - N * ( ";
            if (con0.node->getNameNum() != datum)
                outFile << "x(" << con0.node->getNameNum() << ')';
            if (con1.node->getNameNum() != datum)
                outFile << " - x(" << con1.node->getNameNum() << ')';
            outFile << " ) ) ";
        }
        // TODO：其他情况的处理
        break;

    case VSource:
        if (((con0.node->getNameNum() == nameNum1) && (con0.node->getNameNum() == nameNum2)) || ((con1.node->getNameNum() == nameNum1) && (con1.node->getNameNum() == nameNum2))) {
            JAC[nameNum1][nameNum2] += 1;
            outFile << " + 1";
        } else if (((con0.node->getNameNum() == nameNum1) && (con1.node->getNameNum() == nameNum2)) || ((con1.node->getNameNum() == nameNum1) && (con0.node->getNameNum() == nameNum2))) {
            JAC[nameNum1][nameNum2] -= 1;
            outFile << " - 1";
        }
        break;

    case ISource:
        break;

    case Diode:
        // if (((con0.node->getNum() == nodeNum) && (con0.node->getNameNum() == wrt)) ||
        //     ((con1.node->getNum() == nodeNum) && (con1.node->getNameNum() == wrt)))
        // {
        //     outFile << " (" << name << "IS"
        //             << " * " << name << "N"
        //             << ")*(exp("
        //             << name << "N*(";
        //     if (con0.node->getNameNum() != datum)
        //         outFile << "X(" << con0.node->getNameNum() << ')';
        //     if (con1.node->getNameNum() != datum)
        //         outFile << " - X(" << con1.node->getNameNum() << ')';
        //     outFile << "))) ";
        // }
        // else if (((con0.node->getNum() == nodeNum) && (con1.node->getNameNum() == wrt)) ||
        //          ((con1.node->getNum() == nodeNum) && (con0.node->getNameNum() == wrt)))
        // {
        //     outFile << " ( - " << name << "IS"
        //             << " * " << name << "N"
        //             << ")*(exp("
        //             << name << "N*(";
        //     if (con0.node->getNameNum() != datum)
        //         outFile << "X(" << con0.node->getNameNum() << ')';
        //     if (con1.node->getNameNum() != datum)
        //         outFile << " - X(" << con1.node->getNameNum() << ')';
        //     outFile << "))) ";
        // }
        // else
        //     outFile << " 0 ";
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
        // outFile << " 0 ";
        break;

    case Inductor:
        cerr << "This section is not completed" << endl
             << "PROGRAM ENDED ABNORMALLY!" << endl;
        exit(0);
        break;
    }
    return;
}

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

Node *Node::getNext() {
    return next;
}

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

/*
 * 生成节点分析法的方程
 * @param outFile 输出文件流
 * @param F_x 方程矩阵向量
 * @param X 未知变量向量
 * @param datum 接地节点编号
 * @param lastnode 最大节点编号
 */
void Node::genNodalEquation(ofstream &outFile, vector<double> &F_x, vector<double>& X, int datum, int lastnode) {
    // Connections *conList = getConList();
    // outFile << "F(" << nameNum << ") = ";
    // while(conList != NULL){
    //     conList->comp->genEquation(outFile, F_x, X, datum, lastnode, nameNum);
    //     conList = conList->next;
    // }
    Connections *conList = getConList();
    Boolean isHaveVSource = FALSE;
    int VSourceID = NA; // 记录一下电压源编号
    while(conList != NULL) {
        if(conList->comp->getType() == VSource){
            isHaveVSource = TRUE;
            break;
        }
        conList = conList->next;
    }
    if(isHaveVSource == TRUE) {
        // KVL 方程处理
        // assert：假定除接地节点外，每个节点只能接一个电压源
        // TODO：处理节点接多电压源的情况
        conList = getConList();
        outFile << "F(" << nameNum << ") = ";
        while(conList != NULL) {
            if(conList->comp->getType() == VSource) { // KVL
                conList->comp->genEquation(outFile, F_x, X, datum, lastnode, nameNum);
                outFile << endl;
                VSourceID = conList->comp->getcompNum();
                break;
            }
            conList = conList->next;
        }
        // 补偿方程处理
        conList = getConList();
        // outFile<<endl<<nameNum<<"   "<<conList->comp->getcompNum()<<endl;
        F_x[lastnode + VSourceID] += X[lastnode + VSourceID];
        outFile << "F(" << lastnode + VSourceID << ") = x(" << lastnode + VSourceID << ")";
        while(conList != NULL) {
            if(conList->comp->getType() != VSource) {
                conList->comp->genEquation(outFile, F_x, X, datum, lastnode, nameNum);
            }
            conList = conList->next;
        }
        outFile << endl;
    } else {
        // 没有电压源的情况
        Connections *conList = getConList();
        outFile << "F(" << nameNum << ") = ";
        while(conList != NULL){
            conList->comp->genEquation(outFile, F_x, X, datum, lastnode, nameNum);
            conList = conList->next;
        }
        outFile<<endl;
    }
}

/*
 * 生成节点分析法方程的雅可比矩阵
 * @param outFile 输出文件流
 * @param JAC F_x的雅可比矩阵
 * @param X 未知变量向量
 * @param nodePtr2 二维节点指针
 * @param datum 接地节点编号
 * @param lastnode 最大节点编号
 */
void Node::genNodalJAC(ofstream &outFile, vector<vector<double>> &JAC, vector<double>& X, Node *nodePtr2, int datum, int lastnode) {
    Connections *conList = getConList();
    Boolean isHaveVSource = FALSE;
    int VSourceID = NA;
    while(conList != NULL) {
        if(conList->comp->getType() == VSource){
            isHaveVSource = TRUE;
            VSourceID = conList->comp->getcompNum();
            break;
        }
        conList = conList->next;
    }
    if(isHaveVSource == TRUE) {
        // KVL 方程求导处理
        // assert：假定除接地节点外，每个节点只能接一个电压源
        // TODO：处理节点接多电压源的情况
        conList = getConList();
        outFile << "JAC(" << nameNum << "," << nodePtr2->getNameNum() << ") = ";
        while(conList != NULL) {
            if(conList->comp->getType() == VSource) { // KVL
                conList->comp->genJAC(outFile, JAC, X, nameNum, nodePtr2->getNameNum(), datum, lastnode);
                outFile << endl;
                break;
            }
            conList = conList->next;
        }
        // 补偿方程处理
        conList = getConList();
        double temp = JAC[nameNum][nodePtr2->getNameNum()];
        outFile << "JAC(" << lastnode + VSourceID << "," << nodePtr2->getNameNum() << ") = ";
        while(conList != NULL) {
            if(conList->comp->getType() != VSource) {
                conList->comp->genJAC(outFile, JAC, X, nameNum, nodePtr2->getNameNum(), datum, lastnode);
            }
            conList = conList->next;
        }
        JAC[lastnode + VSourceID][nodePtr2->getNameNum()] += JAC[nameNum][nodePtr2->getNameNum()] - temp;
        JAC[nameNum][nodePtr2->getNameNum()] = temp;
        outFile << endl;
    } else {
        // 没有电压源的情况
        Connections *conList = getConList();
        outFile << "JAC(" << nameNum << "," << nodePtr2->getNameNum() << ") = ";
        while(conList != NULL){
            conList->comp->genJAC(outFile, JAC, X, nameNum, nodePtr2->getNameNum(), datum, lastnode);
            conList = conList->next;
        }
        outFile<<endl;
    }
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
    strcpy(name, nameIn);
    type = typeIn;
    is = isIn;
    bf = bfIn;
    br = brIn;
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

// 获取模型temp参数值
double Model::getTemp() {
    return temp;
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