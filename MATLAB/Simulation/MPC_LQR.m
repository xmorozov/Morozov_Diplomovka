clc;clear;close all

%% State Representation
syms theta0 theta1 d_theta0 d_theta1 d2_theta0  d2_theta1 tau
syms m0 m1 L0 L1 l1 I0 I1 g
alfa = I0+L0^2*m1+l1^2*m1*sin(theta1)^2;
beta = L0*m1*l1*sin(theta1);
gama = L0*m1*l1*cos(theta1);
delta = I1+l1^2*m1;
eps = l1^2*m1*sin(theta1)*cos(theta1);
rho = m1*g*l1*sin(theta1);
sigma = 2*l1^2*m1*sin(theta1)*cos(theta1);

states = [d_theta0;...
         (gama*(eps*d_theta0^2+rho)-delta*(tau+beta*d_theta1^2-sigma*d_theta0*d_theta1))/(gama^2-alfa*delta);...
          d_theta1;...
          (gama*(tau+beta*d_theta1^2-sigma*d_theta0*d_theta1)-alfa*(eps*d_theta1^0+rho))/(gama^2-alfa*delta)];
      
      
A = [diff(states(1,1),theta0), diff(states(1,1),d_theta0), diff(states(1,1),theta1), diff(states(1,1),d_theta1);...
    diff(states(2,1),theta0), diff(states(2,1),d_theta0), diff(states(2,1),theta1), diff(states(2,1),d_theta1);...
    diff(states(3,1),theta0), diff(states(3,1),d_theta0), diff(states(3,1),theta1), diff(states(3,1),d_theta1);...
    diff(states(4,1),theta0), diff(states(4,1),d_theta0), diff(states(4,1),theta1), diff(states(4,1),d_theta1)];

B = [diff(states(1,1),tau);...
    diff(states(2,1),tau);...
    diff(states(3,1),tau);...
    diff(states(4,1),tau)];
   
%% System parametres
m0 = 0.6; % mass of arm
m1 = 0.198; % mass of pendulum
L0 = 0.51; % arm length
L1 = 0.23; % pendulum length
l1 = 0.21; % location of pendulum center mass
I0 = (1/3)*m0*L0^2; % moment of iteria of the arm
I1 = (1/12)*m1*L1^2; % moment of iteria of pendulum
g = 9.81; % gravity

%% Operation points
tau=0;
theta0=0;
theta1=0;
d_theta0=0;
d_theta1=0;
A_up = [0, 1, 0, 0;...
      0, -(2*L0*d_theta0*l1^3*m1^2*cos(theta1)^2*sin(theta1) + 2*d_theta1*l1^2*m1*cos(theta1)*sin(theta1)*(m1*l1^2 + I1))/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2), ((m1*l1^2 + I1)*(L0*d_theta1^2*l1*m1*cos(theta1) - 2*d_theta0*d_theta1*l1^2*m1*cos(theta1)^2 + 2*d_theta0*d_theta1*l1^2*m1*sin(theta1)^2) - L0*l1*m1*cos(theta1)*(g*l1*m1*cos(theta1) + d_theta0^2*l1^2*m1*cos(theta1)^2 - d_theta0^2*l1^2*m1*sin(theta1)^2) + L0*l1*m1*sin(theta1)*(m1*cos(theta1)*sin(theta1)*d_theta0^2*l1^2 + g*m1*sin(theta1)*l1))/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2) - ((2*L0^2*l1^2*m1^2*cos(theta1)*sin(theta1) + 2*l1^2*m1*cos(theta1)*sin(theta1)*(m1*l1^2 + I1))*((m1*l1^2 + I1)*(L0*m1*sin(theta1)*d_theta1^2*l1 - 2*d_theta0*m1*cos(theta1)*sin(theta1)*d_theta1*l1^2 + tau) - L0*l1*m1*cos(theta1)*(m1*cos(theta1)*sin(theta1)*d_theta0^2*l1^2 + g*m1*sin(theta1)*l1)))/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2)^2,      -((m1*l1^2 + I1)*(2*d_theta0*m1*cos(theta1)*sin(theta1)*l1^2 - 2*L0*d_theta1*m1*sin(theta1)*l1))/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2);...
      0, 0, 0, 1;...
      0, (2*L0*d_theta1*l1^3*m1^2*cos(theta1)^2*sin(theta1))/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2), ((l1^2*m1*cos(theta1)^2 - l1^2*m1*sin(theta1)^2 + g*l1*m1*cos(theta1))*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0*l1*m1*cos(theta1)*(L0*d_theta1^2*l1*m1*cos(theta1) - 2*d_theta0*d_theta1*l1^2*m1*cos(theta1)^2 + 2*d_theta0*d_theta1*l1^2*m1*sin(theta1)^2) + L0*l1*m1*sin(theta1)*(L0*m1*sin(theta1)*d_theta1^2*l1 - 2*d_theta0*m1*cos(theta1)*sin(theta1)*d_theta1*l1^2 + tau) + 2*l1^2*m1*cos(theta1)*sin(theta1)*(m1*cos(theta1)*sin(theta1)*l1^2 + g*m1*sin(theta1)*l1))/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2) - ((2*L0^2*l1^2*m1^2*cos(theta1)*sin(theta1) + 2*l1^2*m1*cos(theta1)*sin(theta1)*(m1*l1^2 + I1))*((m1*cos(theta1)*sin(theta1)*l1^2 + g*m1*sin(theta1)*l1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0*l1*m1*cos(theta1)*(L0*m1*sin(theta1)*d_theta1^2*l1 - 2*d_theta0*m1*cos(theta1)*sin(theta1)*d_theta1*l1^2 + tau)))/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2)^2, (L0*l1*m1*cos(theta1)*(2*d_theta0*m1*cos(theta1)*sin(theta1)*l1^2 - 2*L0*d_theta1*m1*sin(theta1)*l1))/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2)];
 
B_up = [0;...
        (m1*l1^2 + I1)/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2);...
        0;...
        -(L0*l1*m1*cos(theta1))/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2)];
theta1=pi; 
A_down = [0, 1, 0, 0;...
      0, -(2*L0*d_theta0*l1^3*m1^2*cos(theta1)^2*sin(theta1) + 2*d_theta1*l1^2*m1*cos(theta1)*sin(theta1)*(m1*l1^2 + I1))/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2), ((m1*l1^2 + I1)*(L0*d_theta1^2*l1*m1*cos(theta1) - 2*d_theta0*d_theta1*l1^2*m1*cos(theta1)^2 + 2*d_theta0*d_theta1*l1^2*m1*sin(theta1)^2) - L0*l1*m1*cos(theta1)*(g*l1*m1*cos(theta1) + d_theta0^2*l1^2*m1*cos(theta1)^2 - d_theta0^2*l1^2*m1*sin(theta1)^2) + L0*l1*m1*sin(theta1)*(m1*cos(theta1)*sin(theta1)*d_theta0^2*l1^2 + g*m1*sin(theta1)*l1))/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2) - ((2*L0^2*l1^2*m1^2*cos(theta1)*sin(theta1) + 2*l1^2*m1*cos(theta1)*sin(theta1)*(m1*l1^2 + I1))*((m1*l1^2 + I1)*(L0*m1*sin(theta1)*d_theta1^2*l1 - 2*d_theta0*m1*cos(theta1)*sin(theta1)*d_theta1*l1^2 + tau) - L0*l1*m1*cos(theta1)*(m1*cos(theta1)*sin(theta1)*d_theta0^2*l1^2 + g*m1*sin(theta1)*l1)))/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2)^2,      -((m1*l1^2 + I1)*(2*d_theta0*m1*cos(theta1)*sin(theta1)*l1^2 - 2*L0*d_theta1*m1*sin(theta1)*l1))/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2);...
      0, 0, 0, 1;...
      0, (2*L0*d_theta1*l1^3*m1^2*cos(theta1)^2*sin(theta1))/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2), ((l1^2*m1*cos(theta1)^2 - l1^2*m1*sin(theta1)^2 + g*l1*m1*cos(theta1))*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0*l1*m1*cos(theta1)*(L0*d_theta1^2*l1*m1*cos(theta1) - 2*d_theta0*d_theta1*l1^2*m1*cos(theta1)^2 + 2*d_theta0*d_theta1*l1^2*m1*sin(theta1)^2) + L0*l1*m1*sin(theta1)*(L0*m1*sin(theta1)*d_theta1^2*l1 - 2*d_theta0*m1*cos(theta1)*sin(theta1)*d_theta1*l1^2 + tau) + 2*l1^2*m1*cos(theta1)*sin(theta1)*(m1*cos(theta1)*sin(theta1)*l1^2 + g*m1*sin(theta1)*l1))/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2) - ((2*L0^2*l1^2*m1^2*cos(theta1)*sin(theta1) + 2*l1^2*m1*cos(theta1)*sin(theta1)*(m1*l1^2 + I1))*((m1*cos(theta1)*sin(theta1)*l1^2 + g*m1*sin(theta1)*l1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0*l1*m1*cos(theta1)*(L0*m1*sin(theta1)*d_theta1^2*l1 - 2*d_theta0*m1*cos(theta1)*sin(theta1)*d_theta1*l1^2 + tau)))/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2)^2, (L0*l1*m1*cos(theta1)*(2*d_theta0*m1*cos(theta1)*sin(theta1)*l1^2 - 2*L0*d_theta1*m1*sin(theta1)*l1))/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2)];
 
B_down = [0;...
        (m1*l1^2 + I1)/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2);...
        0;...
        -(L0*l1*m1*cos(theta1))/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2)];
      
%% System
Ts=0.02;
[nx,nu] = size(B);
C = eye(nx);
ny = size(C,1);
D = zeros(ny,nu);
sys_c = ss(A_up,B_up,C,D);
sys_d = c2d(sys_c,Ts);

% simulations
k = 300;
x0 = [-1;-2;0.5;2];
%x0 = [-0.3151;-2.6642;-2.0901;-11.2293];
%% LQR setup
lqr_x = zeros(nx,k);
lqr_u = zeros(nu,k);
lqr_y = zeros(nx,k);
Q = diag([0.6 0.05 10 0.5]);
R = 1;
K = dlqr(sys_d.A,sys_d.B,Q,R);

%% MPC setup
mpc_x = zeros(nx,k);
mpc_u = zeros(nu,k);
mpc_y = zeros(nx,k);
N = 20;
u = sdpvar(nu,N);
x = sdpvar(nx,N+1);
y = sdpvar(nx,N);
Qx = diag([1.5 0.08 10 0.2]);
Qu = 1;
u_cst = [-5,5];
cst = [];
obj = 0;
for i = 1:N
    obj = obj+x(:,i)'*Qx*x(:,i)+u(:,i)'*Qu*u(:,i);
    cst = cst+[x(:,i+1) == sys_d.A*x(:,i)+sys_d.B*u(:,i)];
    cst = cst+[y(:,i) == sys_d.C*x(:,i)+sys_d.D*u(:,i)];
    cst = cst+[u_cst(1)<=u(:,i)<=u_cst(2)];
end

%% LQR 
lqr_x(:,1) = x0;
for i=1:k
   lqr_u(:,i) = -K*lqr_x(:,i);
   lqr_x(:,i+1) = sys_d.A*lqr_x(:,i)+sys_d.B*lqr_u(:,i);
   lqr_y(:,i) = sys_d.C*lqr_x(:,i);
end

%% MPC 
options = sdpsettings('verbose',0,'cachesolvers',1,'solver','quadprog');
OPT = optimizer(cst,obj,options,x(:,1),u(:,1));
mpc_x(:,1)=x0;
for i = 1:k
    % optimizacia
    mpc_u(:,i) = OPT(mpc_x(:,i));

    % simulacia
    mpc_x(:,i+1) = sys_d.A*mpc_x(:,i)+sys_d.B*mpc_u(:,i);
    mpc_y(:,i) = sys_d.C*mpc_x(:,i)+sys_d.D*mpc_u(:,i);

end

%% Plots
t = linspace(0,k*Ts,k);
figure
txt = {'\theta_0','d\theta_0','\theta_1','d\theta_1'};
for i = 1:nx
    subplot(nx+1,1,i)
    hold on
    plot(t,lqr_y(i,1:end))
    plot(t,mpc_y(i,1:end))
    title(txt{i})
    legend('LQR','MPC')
    xlabel('t [s]')
    hold off
end
subplot(nx+1,1,nx+1)
hold on
plot(t,lqr_u)
plot(t,mpc_u)
hold off
legend('LQR','MPC')
xlabel('t [s]')
title('u')



