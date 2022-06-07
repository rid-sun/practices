#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>

namespace config {
    std::string inFileName;
    std::string outFileName;
    bool help = false;
    bool instrError = false;
    int datum = -1;

    void parse_args(int argc, char** argv);
}  // namespace config

/*
 * 用来解析参数
 * 严格限制命令形式，必须是parser <-f filename> <-o filename> <-h> <-d datum>的形式
 */
void config::parse_args(int argc, char** argv) {
    for (int i = 1; i < argc;i++) {
        if(argv[i][0] == '-') {
            if(argv[i] == std::string("-o")) {
                if(std::string("-") == argv[i + 1])
                    outFileName = ""; 
                else
                    outFileName = argv[i + 1];
                i++;
            }else if(argv[i] == std::string("-f")) {
                if (std::string("-") == argv[i + 1])
                    inFileName = "";
                else
                    inFileName = argv[i + 1];
                i++;
            }else if(argv[i] == std::string("-h")) {
                help = true;
            }else if(argv[i] == std::string("-d")) {
                datum = atoi(argv[i + 1]);
                i++;
            }else {
                instrError = true;
            }
        }else {
            instrError = true;
        }
    }
}