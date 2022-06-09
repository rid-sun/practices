function [F] = equation(X)

% EQUATION Calculates the value of the system of equations for Chua Circuit (Modified Nodal Equations)

% Circuit Parameters

R1 = 10000;
R2 = 4000;
R3 = 4000;
R4 = 5000;
R5 = 30000;
R6 = 500;
R7 = 500;
R8 = 30000;
R9 = 10100;
R10 = 10100;
R11 = 4000;
R12 = 4000;
R13 = 30000;
R14 = 30000;

% Voltage Sources

VCC1 = 12;
VCC2 = 12;
V1 = 10;
V2 = 2;

% BJT Parameters

Q1IS = 1e-9;
Q1BF = 100;
Q1BR = 1;
Q1N = 38.7766;
Q2IS = 1e-9;
Q2BF = 100;
Q2BR = 1;
Q2N = 38.7766;
Q3IS = 1e-9;
Q3BF = 100;
Q3BR = 1;
Q3N = 38.7766;
Q4IS = 1e-9;
Q4BF = 100;
Q4BR = 1;
Q4N = 38.7766;

% Ebers-Moll Transistor Model

Ic1 = Q1IS*(exp(-Q1N*(X(8)-X(1)))-1) - Q1IS/Q1BR*(1+Q1BR)*(exp(-Q1N*(X(4)-X(1)))-1);
Ic2 = Q2IS*(exp(-Q2N*(X(8)-X(6)))-1) - Q2IS/Q2BR*(1+Q2BR)*(exp(-Q2N*(X(10)-X(6)))-1);
Ic3 = Q3IS*(exp(-Q3N*(X(9)-X(1)))-1) - Q3IS/Q3BR*(1+Q3BR)*(exp(-Q3N*(X(5)-X(1)))-1);
Ic4 = Q4IS*(exp(-Q4N*(X(9)-X(7)))-1) - Q4IS/Q4BR*(1+Q4BR)*(exp(-Q4N*(X(12)-X(7)))-1);

Ie1 = - Q1IS/Q1BF*(1+Q1BF)*(exp(-Q1N*(X(8)-X(1)))-1) + Q1IS*(exp(-Q1N*(X(4)-X(1)))-1);
Ie2 = - Q2IS/Q2BF*(1+Q2BF)*(exp(-Q2N*(X(8)-X(6)))-1) + Q2IS*(exp(-Q2N*(X(10)-X(6)))-1);
Ie3 = - Q3IS/Q3BF*(1+Q3BF)*(exp(-Q3N*(X(9)-X(1)))-1) + Q3IS*(exp(-Q3N*(X(5)-X(1)))-1);
Ie4 = - Q4IS/Q4BF*(1+Q4BF)*(exp(-Q4N*(X(9)-X(7)))-1) + Q4IS*(exp(-Q4N*(X(12)-X(7)))-1);

% System of Equations for Chua Circuit (Modified Nodal Equations)

F(1) = (X(1)-X(2)) - V1;

F(2) = (X(2)-X(3))/R1 - X(17);

F(3) = (X(3)-X(2))/R1 + (X(3)-X(10))/R13 - X(18);

F(4) = (X(4)-X(13))/R2 + (X(4)-X(6))/R5 + Ic1;

F(5) = (X(5)-X(14))/R3 + (X(5)-X(7))/R8 + Ic3;

F(6) = (X(6)-X(4))/R5 + X(6)/R9 - Ic2 - Ie2; 

F(7) = (X(7)-X(5))/R8 + X(7)/R10 - Ic4 - Ie4;  

F(8) = X(8)/R6 + Ie1 + Ie2;

F(9) = X(9)/R7 + Ie3 + Ie4;

F(10) = (X(10)-X(13))/R11 + (X(10)-X(3))/R13 + Ic2;

F(11) = (X(11)-X(3)) - V2;

F(12) = (X(12)-X(14))/R12 + (X(12)-X(11))/R14 + Ic4;

F(13) = X(13) - VCC1;

F(14) = X(14) - VCC2;

F(15) = (X(13)-X(4))/R2 + (X(13)-X(10))/R11 + X(15);

F(16) = (X(14)-X(5))/R3 + (X(14)-X(12))/R12 + X(16);

F(17) = X(1)/R4 + X(17) - Ic1 - Ie1 - Ic3 - Ie3;

F(18) = (X(11)-X(12))/R14 + X(18);

F = F';

end