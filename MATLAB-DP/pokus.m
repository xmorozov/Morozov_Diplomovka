%% pokus
clc,clear,close all
%% 

% mpc_y(:,50)
% ans =
%    -0.1604
%     0.7438
%    -0.0282
%    -0.0138
% 
% mpc_y(:,51)
% ans = 
%    -0.1456
%     0.7317
%    -0.0285
%    -0.0122
% 
% mpc_u(:,50) 
% ans =
%    -0.0596
x0 = [-0.1604 0.7438 -0.0282 -0.0138]';
u = -0.0596;
tspan = [0, 0.02];
opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
[t,y] = ode45(@(t,x) odefun(t,x,u),tspan,x0,opts);
















