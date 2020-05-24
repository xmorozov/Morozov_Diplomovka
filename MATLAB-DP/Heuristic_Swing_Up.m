clc;clear;
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

x_states = [d_theta0;...
         (gama*(eps*d_theta0^2+rho)-delta*(tau+beta*d_theta1^2-sigma*d_theta0*d_theta1))/(gama^2-alfa*delta);...
          d_theta1;...
          (gama*(tau+beta*d_theta1^2-sigma*d_theta0*d_theta1)-alfa*(eps*d_theta1^0+rho))/(gama^2-alfa*delta)];
      
      
A = [diff(x_states(1,1),theta0), diff(x_states(1,1),d_theta0), diff(x_states(1,1),theta1), diff(x_states(1,1),d_theta1);...
    diff(x_states(2,1),theta0), diff(x_states(2,1),d_theta0), diff(x_states(2,1),theta1), diff(x_states(2,1),d_theta1);...
    diff(x_states(3,1),theta0), diff(x_states(3,1),d_theta0), diff(x_states(3,1),theta1), diff(x_states(3,1),d_theta1);...
    diff(x_states(4,1),theta0), diff(x_states(4,1),d_theta0), diff(x_states(4,1),theta1), diff(x_states(4,1),d_theta1)];

B = [diff(x_states(1,1),tau);...
    diff(x_states(2,1),tau);...
    diff(x_states(3,1),tau);...
    diff(x_states(4,1),tau)];
   
%% System parametres
m0 = 0.6; % mass of arm
m1 = 0.198; % mass of pendulum
L0 = 0.51; % arm length
L1 = 0.23; % pendulum length
l1 = 0.21; % location of pendulum center mass
I0 = (1/3)*m0*L0^2; % moment of iteria of the arm
I1 = (1/12)*m1*L1^2; % moment of iteria of pendulum
g = 9.81; % gravity
om = sqrt(m1*g*l1/I1);
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
C = [0,0,1,0];
ny = size(C,1);
D = zeros(ny,nu);
%up-right
sys1_c = ss(A_up,B_up,C,D);
sys1_d = c2d(sys1_c,Ts);

A = sys1_d.A;
B = sys1_d.B;
C = sys1_d.C;
D = sys1_d.D;

%downside
sys2_c = ss(A_down,B_down,C,D);
sys2_d = c2d(sys2_c,Ts);

A2 = sys2_d.A;
B2 = sys2_d.B;
C2 = sys2_d.C;
D2 = sys2_d.D;

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


options = sdpsettings('verbose',0,'cachesolvers',1,'solver','gurobi');
OPT = optimizer(cst, obj, options, xx{1}, uu{1});

%% Simulation
tf = 15; % sekundy
kf = tf/Ts;
x0 = [0,0,0.01,-0.01];

x_sim = zeros(nx,kf);
u_sim = zeros(nu,kf);
y_sim = zeros(1,kf);
x_sim(:,1)=x0;

IAE = 0;
swg = 1;
a = (m1*g*l1/2);
om = sqrt(m1*g*l1/I1);
opts = odeset('RelTol',1e-2,'AbsTol',1e-4);

for k = 1:kf    
    if swg == 1 %swing up
        E = a*((x_sim(4,k)/om)^2+cos(x_sim(3,k)-1));
        u_sim(:,k) = 4*E*sign(x_sim(4,k)*cos(x_sim(3,k)));
        
        x_sim(:,k+1) = A2*x_sim(:,k)+B2*u_sim(:,k);
        y_sim(k) = C2*x_sim(:,k)+D2*u_sim(:,k)+pi;
        
        x_sim(3,k) = x_sim(3,k) +pi;
        if y_sim(k)<0.9 
            swg = 0;
            u = u_sim(:,k);
            tspan = [0, 0.02];
            [t,y] = ode45(@(t,x) odefun(t,x,u),tspan,x_sim(:,k),opts);
            x_sim(:,k+1) = y(end,:)';

        end
    
    else %MPC
        t0 = cputime;
        u_sim(:,k) = OPT(x_sim(:,k));
        t_solver = cputime - t0
        x_sim(:,k+1) = A*x_sim(:,k)+B*u_sim(:,k);
        y_sim(:,k) = C*x_sim(:,k)+D*u_sim(:,k);

        u = u_sim(:,k);
        tspan = [0, 0.02];
        [t,y] = ode45(@(t,x) odefun(t,x,u),tspan,x_sim(:,k),opts);
        x_sim(:,k+1) = y(end,:)';
    end
    IAE = IAE + abs(x_sim(3,k));
end
x_sim(:,end) = [];
IAE
%% Plots

path = 'Dswing/';
close all
t = linspace(0,kf*Ts,kf);
w = 9;
l = 21;
txt = {'$\theta_0\;[deg]$ ','$\dot{\theta_0}\;[deg s^{-1}]$ ','$\theta_1\;[deg]$','$\dot{\theta_1}\;[deg s^{-1}]$ '};
txt2 = {'arm','darm','pend','dpend'};
for i = 1:nx
    figure(i)
    %subplot(nx+1,1,i)
    hold on
    plot(t,rad2deg(x_sim(i,1:end)))
    ylabel(txt{i},'interpreter','latex')
    xlabel('t [s]')
    hold off
        f2p(txt2{i}, 'Xlim', [0, tf],  'Ytol', 0.05, 'Xtol', 0,...
        'extension', 'pdf','Path', path, 'dpi', 150, 'papersize', [l, w], 'Xsplit', 10,'Ysplit',7);

end


figure(5)
hold on
plot(t,u_sim)
hold off
ylabel('$\tau\;[V]$','interpreter','latex')
xlabel('t [s]')
f2p('control', 'Xlim', [0, tf], 'Ytol', 0.05, 'Xtol', 0,...
'extension', 'pdf', 'dpi', 150, 'Path', path, 'papersize', [l, w], 'Xsplit', 10,'Ysplit',7);
