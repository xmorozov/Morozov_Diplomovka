function [sys,x0,str,ts] = furuta(t,x,u,flag)

switch flag,
case 0
   [sys,x0,str,ts] = mdlInitializeSizes; % Inicializacia
case 1
   sys = mdlDerivatives(t,x,u); % vypocet derivacii
case 3
   sys = mdlOutputs(t,x,u); % vypocet vystupov
case {2, 4, 9} % nepouzite hodnoty premennej flag
   sys = [];
otherwise
   error(['unhandled flag = ',num2str(flag)]); % chybove hlasenie
end

% do tejto funkcie vkladame vlastne udaje
function [sys,x0,str,ts] = mdlInitializeSizes
sizes = simsizes;
sizes.NumContStates  = 4; % pocet spojitych stavov
sizes.NumDiscStates  = 0; % pocet diskretnych stavov
sizes.NumOutputs     = 4; % pocet vystupov
sizes.NumInputs      = 1; % pocet vstupov
sizes.DirFeedthrough = 0; % = 0 v pripade, ze matrica D is nulova, t.j. ak vystupy su definovane len pomocou stavov, inak =1
sizes.NumSampleTimes = 1; % 1 pre spojite systemy
sys = simsizes(sizes); 
x0 = [-1;-2;0.5;2]; % zaciatocna podmienka pre dif. rovnicu
str = []; % str je prazdna matica
ts = [0 0]; % velkost periody vzorkovania, pre spojite systemy 0 0

% do tejto funkcie vkladame vlastne udaje - rovnice dynamiky

function sys = mdlDerivatives(t,x,u)
m0 = 0.6; % mass of arm
m1 = 0.198; % mass of pendulum
L0 = 0.51; % arm length
L1 = 0.23; % pendulum length
l1 = 0.21; % location of pendulum center mass
I0 = (1/3)*m0*L0^2; % moment of iteria of the arm
I1 = (1/12)*m1*L1^2; % moment of iteria of pendulum
g = 9.81; % gravity
    
theta0 = x(1);
d_theta0 = x(2);
theta1 = x(3);
d_theta1 = x(4);

alfa = I0+L0^2*m1+l1^2*m1*sin(theta1)^2;
beta = L0*m1*l1*sin(theta1);
gama = L0*m1*l1*cos(theta1);
delta = I1+l1^2*m1;
eps = l1^2*m1*sin(theta1)*cos(theta1);
rho = m1*g*l1*sin(theta1);
sigma = 2*l1^2*m1*sin(theta1)*cos(theta1); 

sys(1) = d_theta0;
sys(2) = (gama*(eps*d_theta0^2+rho)-delta*(u+beta*d_theta1^2-sigma*d_theta0*d_theta1))/(gama^2-alfa*delta);
sys(3) = d_theta1;
sys(4) = (gama*(u+beta*d_theta1^2-sigma*d_theta0*d_theta1)-alfa*(eps*d_theta1^0+rho))/(gama^2-alfa*delta);


% do tejto funkcie vkladame vlastne udaje - rovnice vystupu
function sys = mdlOutputs(t,x,u)
sys(1) = x(1);
sys(2) = x(2);
sys(3) = x(3);
sys(4) = x(4);






