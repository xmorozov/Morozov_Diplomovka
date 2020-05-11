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
          (gama*(tau+beta*d_theta1^2-sigma*d_theta0*d_theta1)-alfa*(eps*d_theta0^2+rho))/(gama^2-alfa*delta)];
      
      
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
0,                        -(2*L0*d_theta0*l1^3*m1^2*cos(theta1)^2*sin(theta1) + 2*d_theta1*l1^2*m1*cos(theta1)*sin(theta1)*(m1*l1^2 + I1))/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2),                                                                                                                                                                                   ((m1*l1^2 + I1)*(L0*d_theta1^2*l1*m1*cos(theta1) - 2*d_theta0*d_theta1*l1^2*m1*cos(theta1)^2 + 2*d_theta0*d_theta1*l1^2*m1*sin(theta1)^2) - L0*l1*m1*cos(theta1)*(g*l1*m1*cos(theta1) + d_theta0^2*l1^2*m1*cos(theta1)^2 - d_theta0^2*l1^2*m1*sin(theta1)^2) + L0*l1*m1*sin(theta1)*(m1*cos(theta1)*sin(theta1)*d_theta0^2*l1^2 + g*m1*sin(theta1)*l1))/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2) - ((2*L0^2*l1^2*m1^2*cos(theta1)*sin(theta1) + 2*l1^2*m1*cos(theta1)*sin(theta1)*(m1*l1^2 + I1))*((m1*l1^2 + I1)*(L0*m1*sin(theta1)*d_theta1^2*l1 - 2*d_theta0*m1*cos(theta1)*sin(theta1)*d_theta1*l1^2 + tau) - L0*l1*m1*cos(theta1)*(m1*cos(theta1)*sin(theta1)*d_theta0^2*l1^2 + g*m1*sin(theta1)*l1)))/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2)^2,      -((m1*l1^2 + I1)*(2*d_theta0*m1*cos(theta1)*sin(theta1)*l1^2 - 2*L0*d_theta1*m1*sin(theta1)*l1))/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2);...
      0, 0, 0, 1;...
0, (2*d_theta0*l1^2*m1*cos(theta1)*sin(theta1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) + 2*L0*d_theta1*l1^3*m1^2*cos(theta1)^2*sin(theta1))/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2), ((g*l1*m1*cos(theta1) + d_theta0^2*l1^2*m1*cos(theta1)^2 - d_theta0^2*l1^2*m1*sin(theta1)^2)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0*l1*m1*cos(theta1)*(L0*d_theta1^2*l1*m1*cos(theta1) - 2*d_theta0*d_theta1*l1^2*m1*cos(theta1)^2 + 2*d_theta0*d_theta1*l1^2*m1*sin(theta1)^2) + 2*l1^2*m1*cos(theta1)*sin(theta1)*(m1*cos(theta1)*sin(theta1)*d_theta0^2*l1^2 + g*m1*sin(theta1)*l1) + L0*l1*m1*sin(theta1)*(L0*m1*sin(theta1)*d_theta1^2*l1 - 2*d_theta0*m1*cos(theta1)*sin(theta1)*d_theta1*l1^2 + tau))/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2) - (((m1*cos(theta1)*sin(theta1)*d_theta0^2*l1^2 + g*m1*sin(theta1)*l1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0*l1*m1*cos(theta1)*(L0*m1*sin(theta1)*d_theta1^2*l1 - 2*d_theta0*m1*cos(theta1)*sin(theta1)*d_theta1*l1^2 + tau))*(2*L0^2*l1^2*m1^2*cos(theta1)*sin(theta1) + 2*l1^2*m1*cos(theta1)*sin(theta1)*(m1*l1^2 + I1)))/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2)^2, (L0*l1*m1*cos(theta1)*(2*d_theta0*m1*cos(theta1)*sin(theta1)*l1^2 - 2*L0*d_theta1*m1*sin(theta1)*l1))/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2)];
 
B_up = [0;...
          (m1*l1^2 + I1)/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2);
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

A = sys_d.A;
B = sys_d.B;
C = sys_d.C;
D = sys_d.D;

%% MPC setup

N = 20;
Qx = diag([8 0.08 10 0.2]);
Qu = 1;

% vytvaranie sdpvarov
xx = cell(N+1, 1);
uu = cell(N, 1);
yy = cell(N, 1);
for kf = 1:1:N+1
    xx{kf} = sdpvar(nx, 1, 'full');
    if kf<N+1
        uu{kf} = sdpvar(nu, 1, 'full');
        yy{kf} = sdpvar(ny, 1, 'full');
    end
end

% vytvaranie obj a cst
cst = [];
obj = 0;
u_cst = [-5,5];
for i = 1:N
    obj = obj + xx{i}'*Qx*xx{i} + uu{i}'*Qu*uu{i};
    
    cst = cst + [xx{i+1} == A*xx{i} + B*uu{i}];
    cst = cst + [yy{i} == C*xx{i}+D*uu{i}];
    
    cst = cst + [u_cst(1) <= uu{i} <= u_cst(2)];
end


%% Simulacie
tf = 3; % sekundy
kf = tf/Ts;
x0 = [-1;-2;0.5;-0.7];
%x0 = [-0.3151;-2.6642;-2.0901;-11.2293];

options = sdpsettings('verbose',0,'cachesolvers',1,'solver','gurobi');
OPT = optimizer(cst, obj, options, xx{1}, uu{1});

mpc_x = zeros(nx,kf+1);
mpc_u = zeros(nu,kf);
mpc_y = zeros(nx,kf);
mpc_x(:,1)=x0;

opts = odeset('RelTol',1e-2,'AbsTol',1e-4);

for i = 1:kf
    % optimizacia
    t0 = cputime;
    [mpc_u(:,i), problem] = OPT(mpc_x(:,i));
    t_solver = cputime - t0
    % simulations
%     mpc_x(:,i+1) = A*mpc_x(:,i)+B*mpc_u(:,i);
%     mpc_y(:,i) = C*mpc_x(:,i)+D*mpc_u(:,i);
 
    u = mpc_u(:,i);
    tspan = [0, 0.02];
    [t,y] = ode45(@(t,x) odefun(t,x,u),tspan,mpc_x(:,i),opts);
    mpc_x(:,i+1) = y(end,:)';
end
mpc_x(:,end) = [];

%% Plots
t = linspace(0,kf*Ts,kf);
figure
txt = {'\theta_0','d\theta_0','\theta_1','d\theta_1'};
for i = 1:nx
    subplot(nx+1,1,i)
    hold on
    plot(t,mpc_x(i,1:end))
    title(txt{i})
    legend('MPC')
    xlabel('t [s]')
    hold off
end
subplot(nx+1,1,nx+1)
hold on
plot(t,mpc_u)
hold off
legend('MPC')
xlabel('t [s]')
title('u')



