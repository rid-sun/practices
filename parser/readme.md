# PARSER

v3：支持包含 **电阻**、**电压源**、**电流源**、**BJT** 器件的网表识别，并且迭代求解出结果，可正常work

### 待补充内容
    需要一个处理器件送值后，如何转换单位的工具函数
    对于瞬态分析中的电路要求太苛刻，后续需要拓展程序适配更多电路
### 可能遇到问题
1. 出现“由于找不到python39.dll而无法运行程序
    可以在conda的环境中找到这个dll然后复制到build目录下
2. 出现如下图的提示
![error1](pic/error1.png)
    那么可以在环境变量中进行配置，并且将对应变量的值设置为conda目录，(**必须注意重启电脑生效**)，如下
![solution1](pic/solution1.png)
3. 因为在plot.py中引入了matplotlib这个包，如果运行失败的话，那可能是没有install该包的原因，下载即可
    具体的测试方法可以用PyRun_SimpleString("import cv2")语句来逐条测试一下有无报错
4. 出现如下图的提示
![error2](pic/error2.png)
    那么可以做如下的配置，注意配置完后重启cmd。
![solution2](pic/solution2.png)
5. 注意路径的正确性