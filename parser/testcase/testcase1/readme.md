# testcase_1
## Circuit
<div align=center>
    <img src = ../../pic/circuit1.png width = 50%>
</div>

> [`Netlist1.txt`](Netlist1.txt)为本电路的PSPICE格式网表表示

## Circuit Parameters

$R3=10k\Omega$， $Rc1=2k\Omega$， $Rc2=1k\Omega$， $Re=100\Omega$， $Rc1=2k\Omega$， $Vcc=10V$， $Vin=1.5V$  
$m_ea_f=m_ca_r=Is=-1.0*10^{-16}A$， $a_f=0.99$， $a_r=0.5$， $n=-38.78$ $1/V$  

## Ebers-Moll Transistor Model
$Ie_1=fe_1-a_r\ast fc_1$  
$Ic_1=fc_1-a_f\ast fe_1$  
$Ie_2=fe_2-a_r\ast fc_2$  
$Ic_2=fc_2-a_f\ast fe_2$  
$fe_1=\frac{Is}{a_f}\ast (e^{(n\ast (x_{2}-x_{5}))}-1)$  
$fc_1=\frac{Is}{a_r}\ast (e^{(n\ast (x_1-x_5))}-1)$  
$fe_2=\frac{Is}{a_f}\ast (e^{(n\ast (x_2-x_4))}-1)$  
$fc_2=\frac{Is}{a_r}\ast (e^{(n\ast (x_3-x_4))}-1)$  

## Modified Nodal Equations

$F_1 = \frac{x_{1}-x_{4}}{R3}+\frac{x_{1}-x_{6}}{R c 1}+I c_{1}=0$  
$F_2 = \frac{x_{2}}{R e}+Ie_1+Ie_2=0$  
$F_3 = \frac{x_{3}-x_{6}}{Rc2}+Ic_2=0 $  
$F_4 = \frac{x_{4}-x_{1}}{R3}-Ic_2-Ie_2=0 $  
$F_5 = x_3-Vin=0 $  
$F_6 = x_6-Vcc=0 $  
$F_7 = \frac{x_{6}-x_{1}}{Rc 1}+\frac{x_{6}-x_{3}}{Rc 2}+x_7=0 $  
$F_8 = x_8-Ic_1-Ie_1=0 $  

> 详见文件[`equation1.m`](equation1.m)

## Jacobian Equations

> 详见文件[`jacobian1.m`](jacobian1.m)

## Standard Solution

<table>
    <tr>
        <td align="center"> <img src = ../../pic/std_sol1_fig.png> 同伦法求解图像 </td>
        <td align="center"> <img src = ../../pic/std_sol1.png> 同伦法解 vs. PSPICE </td>
    </tr>
</table>

## Test Method
将标准解，进行些许修改，以此作为牛顿迭代点额初值来进行求解，即可验证本项目代码的正确性。
