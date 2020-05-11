clc;clear;%close all

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
          (gama*(tau+beta*d_theta1^2-sigma*d_theta0*d_theta1)-alfa*(eps*d_theta0^2+rho))/(gama^2-alfa*delta)];
      
      
A = [diff(x_states(1,1),theta0), diff(x_states(1,1),d_theta0), diff(x_states(1,1),theta1), diff(x_states(1,1),d_theta1);...
    diff(x_states(2,1),theta0), diff(x_states(2,1),d_theta0), diff(x_states(2,1),theta1), diff(x_states(2,1),d_theta1);...
    diff(x_states(3,1),theta0), diff(x_states(3,1),d_theta0), diff(x_states(3,1),theta1), diff(x_states(3,1),d_theta1);...
    diff(x_states(4,1),theta0), diff(x_states(4,1),d_theta0), diff(x_states(4,1),theta1), diff(x_states(4,1),d_theta1)]

B = [diff(x_states(1,1),tau);...
    diff(x_states(2,1),tau);...
    diff(x_states(3,1),tau);...
    diff(x_states(4,1),tau)]
   
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
0, (2*d_theta0*l1^2*m1*cos(theta1)*sin(theta1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) + 2*L0*d_theta1*l1^3*m1^2*cos(theta1)^2*sin(theta1))/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2), ((g*l1*m1*cos(theta1) + d_theta0^2*l1^2*m1*cos(theta1)^2 - d_theta0^2*l1^2*m1*sin(theta1)^2)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0*l1*m1*cos(theta1)*(L0*d_theta1^2*l1*m1*cos(theta1) - 2*d_theta0*d_theta1*l1^2*m1*cos(theta1)^2 + 2*d_theta0*d_theta1*l1^2*m1*sin(theta1)^2) + 2*l1^2*m1*cos(theta1)*sin(theta1)*(m1*cos(theta1)*sin(theta1)*d_theta0^2*l1^2 + g*m1*sin(theta1)*l1) + L0*l1*m1*sin(theta1)*(L0*m1*sin(theta1)*d_theta1^2*l1 - 2*d_theta0*m1*cos(theta1)*sin(theta1)*d_theta1*l1^2 + tau))/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2) - (((m1*cos(theta1)*sin(theta1)*d_theta0^2*l1^2 + g*m1*sin(theta1)*l1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0*l1*m1*cos(theta1)*(L0*m1*sin(theta1)*d_theta1^2*l1 - 2*d_theta0*m1*cos(theta1)*sin(theta1)*d_theta1*l1^2 + tau))*(2*L0^2*l1^2*m1^2*cos(theta1)*sin(theta1) + 2*l1^2*m1*cos(theta1)*sin(theta1)*(m1*l1^2 + I1)))/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2)^2, (L0*l1*m1*cos(theta1)*(2*d_theta0*m1*cos(theta1)*sin(theta1)*l1^2 - 2*L0*d_theta1*m1*sin(theta1)*l1))/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2)]
 
B_up = [0;...
          (m1*l1^2 + I1)/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2);
        0;...
 -(L0*l1*m1*cos(theta1))/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2)]
theta1=pi; 
A_down = [0, 1, 0, 0;...
0,                        -(2*L0*d_theta0*l1^3*m1^2*cos(theta1)^2*sin(theta1) + 2*d_theta1*l1^2*m1*cos(theta1)*sin(theta1)*(m1*l1^2 + I1))/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2),                                                                                                                                                                                   ((m1*l1^2 + I1)*(L0*d_theta1^2*l1*m1*cos(theta1) - 2*d_theta0*d_theta1*l1^2*m1*cos(theta1)^2 + 2*d_theta0*d_theta1*l1^2*m1*sin(theta1)^2) - L0*l1*m1*cos(theta1)*(g*l1*m1*cos(theta1) + d_theta0^2*l1^2*m1*cos(theta1)^2 - d_theta0^2*l1^2*m1*sin(theta1)^2) + L0*l1*m1*sin(theta1)*(m1*cos(theta1)*sin(theta1)*d_theta0^2*l1^2 + g*m1*sin(theta1)*l1))/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2) - ((2*L0^2*l1^2*m1^2*cos(theta1)*sin(theta1) + 2*l1^2*m1*cos(theta1)*sin(theta1)*(m1*l1^2 + I1))*((m1*l1^2 + I1)*(L0*m1*sin(theta1)*d_theta1^2*l1 - 2*d_theta0*m1*cos(theta1)*sin(theta1)*d_theta1*l1^2 + tau) - L0*l1*m1*cos(theta1)*(m1*cos(theta1)*sin(theta1)*d_theta0^2*l1^2 + g*m1*sin(theta1)*l1)))/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2)^2,      -((m1*l1^2 + I1)*(2*d_theta0*m1*cos(theta1)*sin(theta1)*l1^2 - 2*L0*d_theta1*m1*sin(theta1)*l1))/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2);...
      0, 0, 0, 1;...
0, (2*d_theta0*l1^2*m1*cos(theta1)*sin(theta1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) + 2*L0*d_theta1*l1^3*m1^2*cos(theta1)^2*sin(theta1))/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2), ((g*l1*m1*cos(theta1) + d_theta0^2*l1^2*m1*cos(theta1)^2 - d_theta0^2*l1^2*m1*sin(theta1)^2)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0*l1*m1*cos(theta1)*(L0*d_theta1^2*l1*m1*cos(theta1) - 2*d_theta0*d_theta1*l1^2*m1*cos(theta1)^2 + 2*d_theta0*d_theta1*l1^2*m1*sin(theta1)^2) + 2*l1^2*m1*cos(theta1)*sin(theta1)*(m1*cos(theta1)*sin(theta1)*d_theta0^2*l1^2 + g*m1*sin(theta1)*l1) + L0*l1*m1*sin(theta1)*(L0*m1*sin(theta1)*d_theta1^2*l1 - 2*d_theta0*m1*cos(theta1)*sin(theta1)*d_theta1*l1^2 + tau))/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2) - (((m1*cos(theta1)*sin(theta1)*d_theta0^2*l1^2 + g*m1*sin(theta1)*l1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0*l1*m1*cos(theta1)*(L0*m1*sin(theta1)*d_theta1^2*l1 - 2*d_theta0*m1*cos(theta1)*sin(theta1)*d_theta1*l1^2 + tau))*(2*L0^2*l1^2*m1^2*cos(theta1)*sin(theta1) + 2*l1^2*m1*cos(theta1)*sin(theta1)*(m1*l1^2 + I1)))/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2)^2, (L0*l1*m1*cos(theta1)*(2*d_theta0*m1*cos(theta1)*sin(theta1)*l1^2 - 2*L0*d_theta1*m1*sin(theta1)*l1))/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2)]
 
B_down = [0;...
          (m1*l1^2 + I1)/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2);
        0;...
 -(L0*l1*m1*cos(theta1))/((m1*l1^2 + I1)*(m1*L0^2 + m1*l1^2*sin(theta1)^2 + I0) - L0^2*l1^2*m1^2*cos(theta1)^2)]
     
%% System
Ts=0.02;
[nx,nu] = size(B);
C = [0,0,1,0];
ny = size(C,1);
D = zeros(ny,nu);

% upright
sys1_c = ss(A_up,B_up,C,D);
sys1_d = c2d(sys1_c,Ts);

A1 = sys1_d.A
B1 = sys1_d.B
C1 = sys1_d.C;
D1 = sys1_d.D;
%downside
sys2_c = ss(A_down,B_down,C,D);
sys2_d = c2d(sys2_c,Ts);

A2 = sys2_d.A;
B2 = sys2_d.B;
C2 = sys2_d.C;
D2 = sys2_d.D;

%% Simulation
tf = 10; % sekundy
kf = tf/Ts;
x0 = [0,0,pi,-0.01];

x_sim = zeros(nx,kf);
u_sim = zeros(nu,kf);
y_sim = zeros(1,kf);

x_sim(:,1)=x0;
om = sqrt(m1*g*l1/I1);

opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
for k = 1:kf    
    % Energy_Shaping
        % Correct
    E = (m1*g*l1/2)*((x_sim(4,k)/om)^2+cos(x_sim(3,k))-1);
    u_sim(:,k) = -5*E*sign(x_sim(4,k)*cos(x_sim(3,k)));
        % Error
    E = (m1*g*l1/2)*((x_sim(4,k)/om)^2+cos(x_sim(3,k)-1));
    u_sim(:,k) = 4*E*sign(x_sim(4,k)*cos(x_sim(3,k)));
%     
    % Exponential
    %u_sim(:,k) = 0.055*abs(x_sim(3,k)^4)*sign(x_sim(4,k)*cos(x_sim(3,k)));
    
    % Simulations
        %Lin model
%     x_sim(:,k+1) = A2*x_sim(:,k)+B2*u_sim(:,k);
%     y_sim(k) = C2*x_sim(:,k)+D2*u_sim(:,k)+pi;

        %Nonlin model
    u = u_sim(:,k);
    tspan = [0, 0.02];
    SOL = ode45(@(t,x) odefun(t,x,u),tspan,x_sim(:,k),opts);
    ODE_VAL = deval(tspan,SOL);
    x_sim(:,k+1) = ODE_VAL(:,end);
    y_sim(k) = x_sim(3,k);
end
x_sim(:,end) = [];

%% Plots Short
t = linspace(0,tf,kf);
% subplot(3,1,1)
% plot(t,rad2deg(x_sim(1,:)))
% hold on
% legend({'MPC'},'Orientation','horizontal')
% grid on
% box on
% xlabel('t [s]')
% ylabel('$\theta_0\,[deg]$','interpreter','latex')
% % axis([0 tf -140 30])
% 
% subplot(3,1,2)
% plot(t,rad2deg(y_sim))
% hold on
% legend({'MPC'},'Orientation','horizontal')
% %title('Pendulum position')
% xlabel('t [s]')
% ylabel('$\theta_1\,[deg]$','interpreter','latex')
% % axis([0 tf -40 300])
% grid on
% 
% subplot(3,1,3)
% hold on
% plot(t,u_sim)
% % plot([t(1),t(end)],[u_cst(1),u_cst(1)],'k--');
% % plot([t(1),t(end)],[u_cst(2),u_cst(2)],'k--');
% %title('Torque')
% xlabel('t [s]')
% ylabel('$\tau\,[V]$','interpreter','latex')
% axis([0 tf -5.3 5.3])
% grid on
% box on
% hold off
% legend({'MPC'},'Orientation','horizontal')

% Plost Full

txt = {'\theta_0','d\theta_0','\theta_1','d\theta_1'};
for i = 1:nx
    subplot(nx+1,1,i)
    hold on
    plot(t,x_sim(i,1:end))
    title(txt{i})
    xlabel('t [s]')
    hold off
end
subplot(nx+1,1,nx+1)
hold on
plot(t,u_sim)
xlabel('t [s]')
title('u')

%ISE_M
