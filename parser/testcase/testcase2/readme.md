# testcase_2
## Circuit

<div align=center>
    <img src = ../../pic/circuit2.png width = 50%>
</div>

> [`Netlist2.txt`](Netlist2.txt)为本电路的PSPICE格式网表表示

## Circuit Parameters

$R1=10k\Omega $， $R2=5k\Omega $， $R3=1.25k\Omega $， $R4=1M\Omega $， $Rc1=1.5k\Omega $， $Rc2=1k\Omega $， $Re=100\Omega $， $Vcc=10V $  
$m_ea_f=m_ca_r=Is=-1.0*10^{-16}A $， $a_f=0.99 $， $a_r=0.5 $， $n=-38.78$ $1/V $

## Ebers-Moll Transistor Model
$Ie_1=fe_1-a_r\ast fc_1 $  
$Ic_1=fc_1-a_f\ast fe_1 $  
$Ie_2=fe_2-a_r\ast fc_2 $  
$Ic_2=fc_2-a_f\ast fe_2 $  
$fe_1=\frac{Is}{a_f}\ast (e^{(n\ast (x_1-x_5))}-1) $  
$fc_1=\frac{Is}{a_r}\ast (e^{(n\ast (x_2-x_5))}-1) $  
$fe_2=\frac{Is}{a_f}\ast (e^{(n\ast (x_1-x_4))}-1) $  
$fc_2=\frac{Is}{a_r}\ast (e^{(n\ast (x_3-x_4))}-1) $  

## Modified Nodal Equations
$F_1 = \frac{x_{1}}{Re}+Ie_1+Ie_2=0 $  
$F_2 = \frac{x_{2}-x_{4}}{R1}+\frac{x_{2}-x_{6}}{Rc1}+Ic_1=0 $  
$F_3 = \frac{x_{3}-x_{6}}{Rc 2}+Ic_2=0 $  
$F_4 = \frac{x_{4}-x_{2}}{R1} + \frac{x_{4}}{R4}-Ic_2-Ie_2=0 $  
$F_5 = \frac{x_{5}-x_{6}}{R2} + \frac{x_{5}}{R3}-Ic_1-Ie_1=0 $  
$F_6 = x_6-Vcc=0 $  
$F_7 = \frac{x_{6}-x_{2}}{Rc 1}+\frac{x_{6}-x_{3}}{Rc2}+\frac{x_{6}-x_{5}}{R2}+x_7=0 $  

> 详见文件[`equation2.m`](equation2.m)

## Jacobian Equations
> 详见文件[`jacobian2.m`](jacobian2.m)

## Standard Solution
<table>
    <tr>
        <td align=center> <img src = ../../pic/std_sol2_fig.png> 同伦法求解图像 </td>
        <td align=center> <img src = ../../pic/std_sol2.png> 同伦法解 vs. PSPICE </td>
    </tr>
</table>

### Test Method
将标准解，进行些许修改，以此作为牛顿迭代点额初值来进行求解，即可验证本项目代码的正确性。
