function dxdt = odefun(t,x,u)
theta0=x(1);
d_theta0=x(2);
theta1=x(3);
d_theta1=x(4);  

alfa = (1/3)*0.6*0.51^2+0.51^2*0.198+0.21^2*0.198*sin(theta1)^2;
beta = 0.51*0.198*0.21*sin(theta1);
gama = 0.51*0.198*0.21*cos(theta1);
delta = (1/12)*0.198*0.23^2+0.21^2*0.198;
eps = 0.21^2*0.198*sin(theta1)*cos(theta1);
rho = 0.198*9.81*0.21*sin(theta1);
sigma = 2*0.21^2*0.198*sin(theta1)*cos(theta1);

dxdt = zeros(4,1);
dxdt(1) = d_theta0;
dxdt(2) =(gama*(eps*d_theta0^2+rho)-delta*(u+beta*d_theta1^2-sigma*d_theta0*d_theta1))/(gama^2-alfa*delta);
dxdt(3) = d_theta1;
dxdt(4) = (gama*(u+beta*d_theta1^2-sigma*d_theta0*d_theta1)-alfa*(eps*d_theta0^2+rho))/(gama^2-alfa*delta);

end