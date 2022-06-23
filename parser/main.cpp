#include <iostream>
#include <sstream>
#include <vector>
#include "parser_v2.h"
#include "config.h"
#include "evalUtils.h"
using namespace std;

// 分析求解相关
string inFileName, outFileName;
Netlist netlist;

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
        // 1. 参数、变量等初始化
        netlist.setDatum(config::datum);
        inFileName = config::inFileName;
        outFileName = config::outFileName;

        // 2. 解析网表文件，存储信息到nodeList、compList、modelList中
        parseNetList(netlist, inFileName, outFileName);
        
        // 3. 求值准备
        evaluation::preparation(netlist, outFileName);

        // 4. 牛顿迭代求解
        // evaluation::newtonRaphson();

        // 4. 同伦法求解
        evaluation::newtonIterHomo();
    }
    return 0;
}