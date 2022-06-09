function [JAC] = jacobian(x)

% JACOBIAN Calculates the value of the Jacobian for given values of lambda
%          and x for Schmitt Trigger Circuit 1 (Modified Nodal Equations)

% Circuit Parameters

R3 = 10000;
Rc1 = 2000;
Rc2 = 1000;
Re = 100;

% Vcc = 10;
% Vin = 1.5;

% BJT Parameters

Is = 1e-16;
af = 0.99;
ar = 0.5;
n = 38.78;

% Ebers-Moll Transistor Model

% fe1 = -Is/af*(exp(-n*(x(2)-x(5)))-1);
% fc1 = -Is/ar*(exp(-n*(x(1)-x(5)))-1);
% fe2 = -Is/af*(exp(-n*(x(2)-x(4)))-1);
% fc2 = -Is/ar*(exp(-n*(x(3)-x(4)))-1);

% Ie1 = fe1 - ar*fc1;
% Ic1 = fc1 - af*fe1;
% Ie2 = fe2 - ar*fc2;
% Ic2 = fc2 - af*fe2;

% F(1) = (x(1)-x(4))/R3 + (x(1)-x(6))/Rc1 + Ic1; 

JAC(1,1) = 1/R3 + 1/Rc1 + Is*n/ar*(exp(-n*(x(1)-x(5))));
JAC(1,2) = 0 + 0 - Is*n*(exp(-n*(x(2)-x(5)))); 
JAC(1,3) = 0 + 0 + 0;
JAC(1,4) = - 1/R3 + 0 + 0; 
JAC(1,5) = 0 + 0 + Is*n*(exp(-n*(x(2)-x(5)))) - Is*n/ar*(exp(-n*(x(1)-x(5))));
JAC(1,6) = 0 - 1/Rc1 + 0;
JAC(1,7) = 0 + 0 + 0;
JAC(1,8) = 0 + 0 + 0;

% F(2) = x(2)/Re + Ie1 + Ie2;

JAC(2,1) = 0 - Is*n*(exp(-n*(x(1)-x(5)))) + 0;
JAC(2,2) = 1/Re + Is*n/af*(exp(-n*(x(2)-x(5)))) + Is*n/af*(exp(-n*(x(2)-x(4)))); 
JAC(2,3) = 0 + 0 - Is*n*(exp(-n*(x(3)-x(4))));
JAC(2,4) = 0 + 0 - Is*n/af*(exp(-n*(x(2)-x(4)))) + Is*n*(exp(-n*(x(3)-x(4)))); 
JAC(2,5) = 0 - Is*n/af*(exp(-n*(x(2)-x(5)))) + Is*n*(exp(-n*(x(1)-x(5)))) + 0;
JAC(2,6) = 0 + 0 + 0;
JAC(2,7) = 0 + 0 + 0;
JAC(2,8) = 0 + 0 + 0;

% F(3) = (x(3)-x(6))/Rc2 + Ic2;

JAC(3,1) = 0 + 0;
JAC(3,2) = 0 - Is*n*(exp(-n*(x(2)-x(4)))); 
JAC(3,3) = 1/Rc2 + Is*n/ar*(exp(-n*(x(3)-x(4))));
JAC(3,4) = 0 + Is*n*(exp(-n*(x(2)-x(4)))) - Is*n/ar*(exp(-n*(x(3)-x(4)))); 
JAC(3,5) = 0 + 0;
JAC(3,6) = - 1/Rc2 + 0;
JAC(3,7) = 0 + 0 + 0;
JAC(3,8) = 0 + 0 + 0;

% F(4) = (x(4)-x(1))/R3 - Ic2 - Ie2;

JAC(4,1) = -1/R3 + 0 + 0;
JAC(4,2) = 0 - Is*n/af*(exp(-n*(x(2)-x(4)))) + Is*n*(exp(-n*(x(2)-x(4)))); 
JAC(4,3) = 0 - Is*n/ar*(exp(-n*(x(3)-x(4)))) + Is*n*(exp(-n*(x(3)-x(4))));
JAC(4,4) = 1/R3 + Is*n/af*(exp(-n*(x(2)-x(4)))) - Is*n*(exp(-n*(x(3)-x(4)))) - Is*n*(exp(-n*(x(2)-x(4)))) + Is*n/ar*(exp(-n*(x(3)-x(4)))); 
JAC(4,5) = 0 + 0 + 0 + 0;
JAC(4,6) = 0 + 0 + 0 + 0;
JAC(4,7) = 0 + 0 + 0;
JAC(4,8) = 0 + 0 + 0;

% F(5) = x(5) - Vin;

JAC(5,1) = 0 + 0;
JAC(5,2) = 0 + 0; 
JAC(5,3) = 0 + 0;
JAC(5,4) = 0 + 0; 
JAC(5,5) = 1 + 0;
JAC(5,6) = 0 + 0;
JAC(5,7) = 0 + 0;
JAC(5,8) = 0 + 0;

% F(6) = x(6) - Vcc;

JAC(6,1) = 0 + 0;
JAC(6,2) = 0 + 0; 
JAC(6,3) = 0 + 0;
JAC(6,4) = 0 + 0; 
JAC(6,5) = 0 + 0;
JAC(6,6) = 1 + 0;
JAC(5,7) = 0 + 0;
JAC(6,8) = 0 + 0;

% F(7) = (x(6)-x(1))/Rc1 + (x(6)-x(3))/Rc2 + x(7);

JAC(7,1) = - 1/Rc1 + 0 + 0;
JAC(7,2) = 0 + 0 + 0;
JAC(7,3) = 0 - 1/Rc2 + 0;
JAC(7,4) = 0 + 0 + 0;
JAC(7,5) = 0 + 0 + 0;
JAC(7,6) = 1/Rc1 + 1/Rc2 + 0;
JAC(7,7) = 0 + 0 + 1;
JAC(7,8) = 0 + 0 + 0;

% F(8) = x(8) - Ic1 - Ie1;

JAC(8,1) = 0 + Is*n*(exp(-n*(x(1)-x(5)))) - Is*n/ar*(exp(-n*(x(1)-x(5))));
JAC(8,2) = 0 - Is*n/af*(exp(-n*(x(2)-x(5)))) + Is*n*(exp(-n*(x(2)-x(5))));
JAC(8,3) = 0 + 0 + 0;
JAC(8,4) = 0 + 0 + 0;
JAC(8,5) = 0 + Is*n/af*(exp(-n*(x(2)-x(5)))) - Is*n*(exp(-n*(x(1)-x(5)))) - Is*n*(exp(-n*(x(2)-x(5)))) + Is*n/ar*(exp(-n*(x(1)-x(5))));
JAC(8,6) = 0 + 0 + 0;
JAC(8,7) = 0 + 0 + 0;
JAC(8,8) = 1 + 0 + 0;

% F = F';

end