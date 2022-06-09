function [F] = equation(x)

% EQUATION Calculates the value of the system of equations for Schmitt 
%          Trigger Circuit 2 (Modified Nodal Equations)

% Circuit Parameters

R1 = 10000;
R2 = 5000;
R3 = 1250;
R4 = 1000000;
Rc1 = 1500;
Rc2 = 1000;
Re = 100;

Vcc = 10;

% BJT Parameters

Is = 1e-16;
af = 0.99;
ar = 0.5;
n = 38.78;

% Ebers-Moll Transistor Model

fe1 = -Is/af*(exp(-n*(x(1)-x(5)))-1);
fc1 = -Is/ar*(exp(-n*(x(2)-x(5)))-1);
fe2 = -Is/af*(exp(-n*(x(1)-x(4)))-1);
fc2 = -Is/ar*(exp(-n*(x(3)-x(4)))-1);

Ie1 = fe1 - ar*fc1;
Ic1 = fc1 - af*fe1;
Ie2 = fe2 - ar*fc2;
Ic2 = fc2 - af*fe2;

% System of Equations for Schmitt Trigger Circuit 2 (Modified Nodal Equations)

F(1) = x(1)/Re + Ie1 + Ie2;

F(2) = (x(2)-x(4))/R1 + (x(2)-x(6))/Rc1 + Ic1;

F(3) = (x(3)-x(6))/Rc2 + Ic2;

F(4) = (x(4)-x(2))/R1 + x(4)/R4 - Ic2 - Ie2;

F(5) = (x(5)-x(6))/R2 + x(5)/R3 - Ic1 - Ie1;

F(6) = x(6) - Vcc;

F(7) = (x(6)-x(2))/Rc1 + (x(6)-x(3))/Rc2 + (x(6)-x(5))/R2 + x(7);

F = F';

end