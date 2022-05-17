#include<bits/stdc++.h>
using namespace std;

map<double, double> vb, va;//记录多项式中的每一项，格式为 指数 -> 系数

void resolve(string& line){//解析表达式
    double c, e;
    while(!line.empty()){
        char ch[256], op[2];
        int pos = 0, isx = 0, posx;
        sscanf(line.c_str(), "%[+-]%[^+-]", op, ch);//获取一项
        // for (int i = 0; ch[i] != 0;i++)
        //     cout << ch[i];
        // cout << endl;
        while (ch[pos] != 0){//检查该项是否为常数项
            if (ch[pos] == 'x'){
                isx = 1;
                posx = pos;
            }
            pos++;
        }
        if(isx != 0){//处理该项含变量的情况
            if(pos == 1){
                c = e = 1;
            }else if(posx == 0){
                c = 1;
                sscanf(ch, "x%lf", &e);
            }else if(posx == pos - 1){
                e = 1;
                sscanf(ch, "%lfx", &c);
            }else{
                sscanf(ch, "%lfx%lf", &c, &e);
            }
        }else{//常数项情况
            sscanf(ch, "%lf", &c);
            e = 0;
        }
        if(op[0] == '-')//系数符号判断
            vb[e] -= c;
        else
            vb[e] += c;
        line.erase(0, pos + 1);//删除已判定项
    }
    for(auto i:vb){
        // cout << i.second << " " << i.first << endl;
        if(i.first == 0) continue;//常数项求导后为0，不需考虑
        va[i.first - 1] = i.first * i.second;
    }
    // cout << va.size() << " " << vb.size() << endl;
}

double calculate(double x0){//求解下一项
    double x1 = 0, fx0 = 0, fx0_ = 0;
    for(auto i:vb){
        fx0 += i.second * pow(x0, i.first);
    }
    for(auto i:va){
        fx0_ += i.second * pow(x0, i.first);
    }
    x1 = x0 - fx0 / fx0_;
    return x1;
}

/*
	1.只适合一元多次方程的零点求解
	2.如f(x) = x - 4 + 2x^(1/3)的输入格式为：+x-4+2x0.333333
	(3.)首项补充符号为处理方便
*/
int main(){
    // freopen("in.txt", "r", stdin);
    double e, x0, x1;
    string expression;
    cout << "输入函数表达式f(x): " << endl;
    cin >> expression;
    cout << "输入求值精度epson及给定初值x0: " << endl;
    cin >> e >> x0;
    resolve(expression);
    x1 = calculate(x0);
    while(abs(x1 - x0) > e){//牛顿迭代
        x0 = x1;
        x1 = calculate(x0);
    }
    cout << "近似值为： " << x1 << endl;
    return 0;
}