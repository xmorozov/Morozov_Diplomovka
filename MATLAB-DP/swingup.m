clear
clc
close all
%% State Representation
syms theta0 theta1 d_theta0 d_theta1 tau

m0 = 0.6; % mass of arm
m1 = 0.198; % mass of pendulum
L0 = 0.51; % arm length
L1 = 0.23; % pendulum length
l1 = 0.21; % location of pendulum center mass
I0 = (1/3)*m0*L0^2; % moment of iteria of the arm
I1 = (1/12)*m1*L1^2; % moment of iteria of pendulum
g = 9.81; % gravity
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
   

%% Operation points
A_up = double(subs(A, [theta0, theta1, d_theta0, d_theta1, tau], [0 0 0 0 0]));
B_up = double(subs(B, [theta0 theta1 d_theta0 d_theta1 tau], [0 0 0 0 0]));
A_down = double(subs(A, [theta0, theta1, d_theta0, d_theta1, tau], [0 pi 0 0 0]));
B_down = double(subs(B, [theta0 theta1 d_theta0 d_theta1 tau], [0 pi 0 0 0]));

%% Swings
tf = 4;
Ts = 0.02;
kf = tf/Ts;
[nx,nu] = size(B);
C = [0,0,1,0];
ny = size(C,1);
D = zeros(ny,nu);
sys_c = ss(A_down,B_down,C,D);
sys_d = c2d(sys_c,Ts);
x_swing = zeros(nx,kf+1);
u_swing = zeros(1,kf);
y_swing = zeros(1,kf);
x0 = [0;0;0.01;0];
x_swing(:,1) = x0;
a = (m1*g*l1/2);
om = sqrt(m1*g*l1/I1);
opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
for k = 1:kf
    E = a*((x_swing(4,k)/om)^2+cos(x_swing(3,k)-1));
    u_swing(:,k) = 5*E*sign(x_swing(4,k)*cos(x_swing(3,k)));
    %u_swing(:,k) = 0.055*abs(x_swing(3,k)^4)*sign(x_swing(4,k)*cos(x_swing(3,k)));
    x_swing(:,k+1) = sys_d.A*x_swing(:,k)+sys_d.B*u_swing(:,k);
    y_swing(:,k) = sys_d.C*x_swing(:,k) + pi;
    x_swing(3,k)=x_swing(3,k)+pi;
    
%     u = u_swing(:,k);
%     tspan = [0, 0.02];
%     SOL = ode45(@(t,x) odefun(t,x,u),tspan,x_swing(:,k),opts);
%     ODE_VAL = deval(tspan,SOL);
%     x_swing(:,k+1) = ODE_VAL(:,end);
%     y_swing(k) = x_swing(3,k);
end
x_swing(:,end)=[];
%% Plots
% t = linspace(0,kf*Ts,kf);
% figure(1)
% subplot(2,1,1)
% hold on
% plot(t,y_swing)
% % plot([0 kf*Ts],[-5 -5],"r--")
% % plot([0 kf*Ts],[5 5],"r--")
% axis([0 Tsim -5 10])
% hold off
% subplot(2,1,2)
% plot(t,u_swing)

%% Plots
path = 'swings/';
close all
t = linspace(0,kf*Ts,kf);
w = 9;
l = 42;
txt = {'$\theta_0\;\rm{[deg]}$ ','$\dot{\theta_0}\;[\rm{deg\:s^{-1}}]$ ','$\theta_1\;\rm{[deg]}$','$\dot{\theta_1}\;[\rm{deg\: s^{-1}}]$ '};
txt2 = {'arm','darm','pend','dpend'};
txt3 = {'Position of the Arm','darm','Position of the Pendulum','dpend'};
for i = 1:nx
    figure
    %subplot(nx+1,1,i)
    hold on
    plot(t,rad2deg(x_swing(i,1:end)))
    ylabel(txt{i},'interpreter','latex','FontSize',25)
    title(txt3{i},'FontSize',25)
    xlabel('$\rm{t [s]}$','interpreter','latex','FontSize',25)
    hold off
        f2p(txt2{i}, 'Xlim', [0, tf],  'Ytol', 0.05, 'Xtol', 0,...
        'extension', 'pdf','Path', path, 'dpi', 150, 'papersize', [l, w], 'Xsplit', 8,'Ysplit',4);

end
%subplot(nx+1,1,nx+1)
hold on
figure
plot(t,u_swing)
hold off
ylabel('$\tau\;\rm{[N\:m]}$','interpreter','latex','FontSize',25,'position',[-0.13 5.7220e-06 -1])
xlabel('$\rm{t [s]}$','interpreter','latex','FontSize',25)
title('Control Input','FontSize',25)
f2p('control', 'Xlim', [0, tf], 'Ytol', 0.05, 'Xtol', 0,...
'extension', 'pdf', 'dpi', 150, 'Path', path, 'papersize', [l, w], 'Xsplit', 8,'Ysplit',6);







