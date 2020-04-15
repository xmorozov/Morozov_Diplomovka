function x_plus = x_next(x_actual, u_actual)

    m0 = 0.6; % mass of arm
    m1 = 0.198; % mass of pendulum
    L0 = 0.51; % arm length
    L1 = 0.23; % pendulum length
    l1 = 0.21; % location of pendulum center mass
    I0 = (1/3)*m0*L0^2; % moment of iteria of the arm
    I1 = (1/12)*m1*L1^2; % moment of iteria of pendulum
    g = 9.81; % gravity

    theta0 = x_actual(1);
    d_theta0 = x_actual(2);
    theta1 = x_actual(3);
    d_theta1 = x_actual(4);
    u = u_actual;
    
    alfa = I0+L0^2*m1+l1^2*m1*sin(theta1)^2;
    beta = L0*m1*l1*sin(theta1);
    gama = L0*m1*l1*cos(theta1);
    delta = I1+l1^2*m1;
    eps = l1^2*m1*sin(theta1)*cos(theta1);
    rho = m1*g*l1*sin(theta1);
    sigma = 2*l1^2*m1*sin(theta1)*cos(theta1);
    
    x_plus=[d_theta0;...
      (gama*(eps*d_theta0^2+rho)-delta*(u+beta*d_theta1^2-sigma*d_theta0*d_theta1))/(gama^2-alfa*delta);...
       d_theta1;...
      (gama*(u+beta*d_theta1^2-sigma*d_theta0*d_theta1)-alfa*(eps*d_theta1^0+rho))/(gama^2-alfa*delta)];


end