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
sys1_d = ss(A_up,B_up,C,D);
sys1_d = c2d(sys1_d,Ts);
%downside
sys2_c = ss(A_down,B_down,C,D);
sys2_d = c2d(sys2_c,Ts);

%% LQR setup
Q = diag([0.6 0.05 10 0.5]);
R = 1;
K = dlqr(sys1_d.A,sys1_d.B,Q,R);

%% MPC setup
N = 20;
u = sdpvar(nu,N);
x = sdpvar(nx,N+1);
y = sdpvar(nx,N);
Qx = diag([5 0.08 10 0.2]);
Qu = 1;
u_cst = [-5,5];
cst = [];
obj = 0;
ISE_L = 0;
ISE_M = 0;
for i = 1:N
    obj = obj+x(:,i)'*Qx*x(:,i)+u(:,i)'*Qu*u(:,i);
    cst = cst+[x(:,i+1) == sys1_d.A*x(:,i)+sys1_d.B*u(:,i)];
    cst = cst+[y(:,i) == sys1_d.C*x(:,i)+sys1_d.D*u(:,i)];
    cst = cst+[u_cst(1)<=u(:,i)<=u_cst(2)];
end
options = sdpsettings('verbose',0,'cachesolvers',1,'solver','quadprog');
OPT = optimizer(cst,obj,options,x(:,1),u(:,1));

%% Simulation
kf = 300;
x0 = [0,0,0.01,-0.01];

lqr_x = zeros(nx,kf);
lqr_u = zeros(nu,kf);
lqr_y = zeros(1,kf);

x_sim = zeros(nx,kf);
u_sim = zeros(nu,kf);
y_sim = zeros(1,kf);
x_sim(:,1)=x0;
lqr_x(:,1)=x0; %LQR
swg = 1;
for k = 1:kf
    
    if swg == 1 %swing up
        E = (m1*g*l1/2)*((x_sim(4,k)/om)^2+cos(x_sim(3,k)-1));
        u_sim(:,k) = 4*E*sign(x_sim(4,k)*cos(x_sim(3,k)));
        x_sim(:,k+1) = sys2_d.A*x_sim(:,k)+sys2_d.B*u_sim(:,k);
        y_sim(k) = sys2_d.C*x_sim(:,k)+sys2_d.D*u_sim(:,k)+pi;
        
        %LQR
        lqr_u(k) = u_sim(k);
        lqr_x(:,k+1) = x_sim(:,k+1);
        lqr_y(k) = y_sim(k);
        
        if y_sim(k)<0.9 
           swg = 0;
           x_sim(3,k+1)=x_sim(3,k+1)+pi;
           lqr_x(3,k+1)=lqr_x(3,k+1)+pi;
        end
    
    else %MPC
        u_sim(:,k) = OPT(x_sim(:,k));
        x_sim(:,k+1) = sys1_d.A*x_sim(:,k)+sys1_d.B*u_sim(:,k);
        y_sim(:,k) = sys1_d.C*x_sim(:,k)+sys1_d.D*u_sim(:,k);
        
        %LQR
        lqr_u(:,k) = -K*lqr_x(:,k);
        lqr_x(:,k+1) = sys1_d.A*lqr_x(:,k)+sys1_d.B*lqr_u(:,k);
        lqr_y(:,k) = sys1_d.C*lqr_x(:,k);
        
        %ISE
        ISE_L = ISE_L+(lqr_y(:,k))^2;
        ISE_M = ISE_M+(y_sim(:,k))^2;
    end

end
x_sim(:,end) = [];
lqr_x(:,end) = [];
%% Plots
t = linspace(0,kf*Ts,kf);
subplot(3,1,1)
plot(t,rad2deg(x_sim(1,:)))
hold on
plot(t,rad2deg(lqr_x(1,:)))
legend({'MPC','LQR'},'Orientation','horizontal')
grid on
box on
xlabel('t [s]')
ylabel('$\theta_0\,[deg]$','interpreter','latex')
axis([0 6 -140 30])

subplot(3,1,2)
plot(t,rad2deg(y_sim))
hold on
plot(t,rad2deg(lqr_y))
legend({'MPC','LQR'},'Orientation','horizontal')
%title('Pendulum position')
xlabel('t [s]')
ylabel('$\theta_1\,[deg]$','interpreter','latex')
axis([0 6 -40 300])
grid on
hold off
subplot(3,1,3)
hold on
plot(t,u_sim,t,lqr_u)
plot([t(1),t(end)],[u_cst(1),u_cst(1)],'k--');
plot([t(1),t(end)],[u_cst(2),u_cst(2)],'k--');
%title('Torque')
xlabel('t [s]')
ylabel('$\tau\,[V]$','interpreter','latex')
axis([0 6 -5.3 5.3])
grid on
box on
hold off
legend({'MPC','LQR'},'Orientation','horizontal')
% figure
% txt = {'\theta_0','d\theta_0','\theta_1','d\theta_1'};
% for i = 1:nx
%     subplot(nx+1,1,i)
%     hold on
%     plot(t,y_sim(i,1:end))
%     title(txt{i})
%     legend('MPC')
%     xlabel('t [s]')
%     hold off
% end
% subplot(nx+1,1,nx+1)
% hold on
% plot(t,u_sim)
% hold off
% legend('MPC')
% xlabel('t [s]')
% title('u')
ISE_L
ISE_M
