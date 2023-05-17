# testcase_1
### Circuit
<center> <img src = ../../pic/circuit1.png width = 50%> </center>

> `Netlist1.txt`为本电路的PSPICE格式网表表示
### Circuit Parameters
$ R3=10k\Omega $，$ Rc1=2k\Omega $， $ Rc2=1k\Omega $，$ Re=100\Omega $，$ Rc1=2k\Omega $，$ Vcc=10V $，$ Vin=1.5V $
$ m_ea_f=m_ca_r=Is=-1.0*10^{-16}A $，$ a_f=0.99 $，$ a_r=0.5 $，$ n=-38.78 $ $ 1/V $
### Ebers-Moll Transistor Model
$ Ie_1=fe_1-a_r*fc_1 $
$ Ic_1=fc_1-a_f*fe_1 $
$ Ie_2=fe_2-a_r*fc_2 $
$ Ic_2=fc_2-a_f*fe_2 $
$ fe_1=\frac{Is}{a_f}*(e^{(n*(x_2-x_5))}-1) $
$ fc_1=\frac{Is}{a_r}*(e^{(n*(x_1-x_5))}-1) $
$ fe_2=\frac{Is}{a_f}*(e^{(n*(x_2-x_4))}-1) $
$ fc_2=\frac{Is}{a_r}*(e^{(n*(x_3-x_4))}-1) $
### Modified Nodal Equations
$ F_1 = \frac{x_{1}-x_{4}}{R 3}+\frac{x_{1}-x_{6}}{R c 1}+I c_{1}=0 $
$ F_2 = \frac{x_{2}}{R e}+Ie_1+Ie_2=0 $
$ F_3 = \frac{x_{3}-x_{6}}{Rc 2}+Ic_2=0 $
$ F_4 = \frac{x_{4}-x_{1}}{R 3}-Ic_2-Ie_2=0 $
$ F_5 = x_3-Vin=0 $
$ F_6 = x_6-Vcc=0 $
$ F_7 = \frac{x_{6}-x_{1}}{Rc 1}+\frac{x_{6}-x_{3}}{Rc 2}+x_7=0 $
$ F_8 = x_8-Ic_1-Ie_1=0 $
> 详见文件`equation1.m`
### Jacobian Equations
> 详见文件`jacobian1.m`
### Standard Solution
<table>
    <tr>
        <td> <center> <img src = ../../pic/std_sol1_fig.png> 同伦法求解图像 </center> </td>
        <td> <center> <img src = ../../pic/std_sol1.png>
        同伦法解 vs. PSPICE</center> </td>
    </tr>
</table>
### Test Method
将标准解，进行些许修改，以此作为牛顿迭代点额初值来进行求解，即可验证本项目代码的正确性。
