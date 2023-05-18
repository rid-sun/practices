# testcase_3
## Circuit
<div align=center>
    <img src = ../../pic/circuit3.png width = 50%>
</div>

> [`Netlist3.txt`](Netlist3.txt)为本电路的PSPICE格式网表表示

## Circuit Parameters
$R1=10k\Omega $， $R2=R3=4k\Omega $， $R4=5k\Omega $， $R5=R8=30k\Omega $， $R6=R7=0.5k\Omega $， $R9=R10=10.1k\Omega $， $R11=R12=4k\Omega $， $R13=R14=30k\Omega $， $V1=10V $， $V2=2V $， $Vcc1=12V $， $Vcc2=12V $  
$m_ea_f=m_ca_r=Is=-1.0*10^{-9}A $， $a_f=0.99 $， $a_r=0.5 $， $n=-38.7766$ $1/V $

## Ebers-Moll Transistor Model
$Ie_1=fe_1-a_r\ast fc_1 $  
$Ic_1=fc_1-a_f\ast fe_1 $  
$Ie_2=fe_2-a_r\ast fc_2 $  
$Ic_2=fc_2-a_f\ast fe_2 $  
$Ie_3=fe_3-a_r\ast fc_3 $  
$Ic_3=fc_3-a_f\ast fe_3 $  
$Ie_4=fe_4-a_r\ast fc_4 $  
$Ic_4=fc_4-a_f\ast fe_4 $  
$fe_1=\frac{Is}{a_f}\ast (e^{(n\ast (x_8-x_1))}-1) $  
$fc_1=\frac{Is}{a_r}\ast (e^{(n\ast (x_4-x_1))}-1) $  
$fe_2=\frac{Is}{a_f}\ast (e^{(n\ast (x_8-x_6))}-1) $  
$fc_2=\frac{Is}{a_r}\ast (e^{(n\ast (x_{10}-x_6))}-1) $  
$fe_3=\frac{Is}{a_f}\ast (e^{(n\ast (x_9-x_1))}-1) $  
$fc_3=\frac{Is}{a_r}\ast (e^{(n\ast (x_5-x_1))}-1) $  
$fe_4=\frac{Is}{a_f}\ast (e^{(n\ast (x_9-x_7))}-1) $  
$fc_4=\frac{Is}{a_r}\ast (e^{(n\ast (x_{12}-x_7))}-1) $  

## Modified Nodal Equations
$F_1 = x_1 - x_2-V1=0 $  
$F_2 = \frac{x_{2}-x_{3}}{R1}-x_{17}=0 $  
$F_3 = \frac{x_{3}-x_{2}}{R1}+\frac{x_{3}-x_{10}}{R13}-x_{18}=0 $  
$F_4 = \frac{x_{4}-x_{13}}{R2} + \frac{x_{4}-x_{6}}{R5}+Ic_1=0 $  
$F_5 = \frac{x_{5}-x_{14}}{R3} + \frac{x_{5}- x_{7}}{R8}+Ic_3=0 $  
$F_6 = \frac{x_{6}-x_{4}}{R5}+\frac{x_{6}}{R9}-Ic_2-Ie_2=0 $  
$F_7 = \frac{x_{7}-x_{5}}{R8}+\frac{x_{7}}{R10}-Ic_4-Ie_4=0 $  
$F_8 = \frac{x_{8}}{R6}+Ie_1+Ie_2=0 $  
$F_9 = \frac{x_{9}}{R7}+Ie_3+Ie_4=0 $  
$F_{10}=\frac{x_{10}-x_{13}}{R11}+\frac{x_{10}-x_{3}}{R13}+Ic_2=0 $  
$F_{11}=x_{11}-x_{3}-V2=0 $  
$F_{12}=\frac{x_{12}-x_{14}}{R12}+\frac{x_{12}-x_{11}}{R14}+Ic_4=0 $  
$F_{13}=x_{13}-Vcc1=0 $  
$F_{14}=x_{14}-Vcc2=0 $  
$F_{15}=\frac{x_{13}-x_{4}}{R2}+\frac{x_{13}-x_{10}}{R11}+x_{15}=0 $  
$F_{16}=\frac{x_{14}-x_{5}}{R13}+\frac{x_{14}-x_{12}}{R12}+x_{16}=0 $  
$F_{17}=\frac{x_{1}}{R4}+x_{17}-Ic_1-Ie_1-Ic_3-Ie_3=0 $  
$F_{18}=\frac{x_{11}-x_{12}}{R14}+x_{18}=0 $  

> 详见文件[`equation3.m`](equation3.m)  

## Jacobian Equations
> 详见文件[`jacobian3.m`](jacobian3)

## Standard Solution
<div align=center>
    <img src = ../../pic/std_sol3_fig.png width = 50%>
    <p> 同伦法求解图像 </p>
</div>

<table>
    <tr>
        <td align=center><img src = ../../pic/std_sol3_mat.png> 同伦法解 </td>
        <td align=center><img src = ../../pic/std_sol3_ps.png>PSPICE </td>
    </tr>
</table>

## Test Method
将标准解，进行些许修改，以此作为牛顿迭代点额初值来进行求解，即可验证本项目代码的正确性。
