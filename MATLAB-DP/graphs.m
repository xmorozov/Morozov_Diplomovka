clear
clc
%% data setup
load mpc_swup_sim.mat
xx(:,end)=[];
Ts = 0.02;
t = linspace(0,length(xx)*Ts,length(xx));
path = '../';
set(0,'defaulttextinterpreter','latex')
w = 9;
l = 21;
%% plot arm position
figure
hold on
plot([0 20],[0 0],'k-.','LineWidth',1)
plot(t,xx(1,:),'color',[0 0.45 0.74])
xlabel('t')
ylabel('$\theta_0$')
set(get(gca,'ylabel'),'rotation',0)
title('\textbf{Arm position}')
legend('zero','arm position','Orientation','horizontal')
f= get(gcf);
f2p('arm_p', 'Xlim', [0, 20], 'Ylim', [-6.5 1], 'Ytol', 0.05, 'Xtol', 0,...
        'extension', 'pdf', 'Path', path, 'dpi', 150, 'papersize', [l, w], 'Xsplit', 4,'Ysplit',3);
hold off  
%% plot arm velocity
figure
hold on
plot([0 20],[0 0],'k-.','LineWidth',1)
plot(t,xx(2,:),'color',[0 0.45 0.74])
xlabel('t')
ylabel('$\dot{\theta}_0$')
set(get(gca,'ylabel'),'rotation',0)
title('\textbf{Arm angular velocity}')
legend('zero','arm velocity','Orientation','horizontal')
f= get(gcf);
f2p('arm_v', 'Xlim', [0, 20], 'Ylim', [-10 6.5], 'Ytol', 0.05, 'Xtol', 0,...
        'extension', 'pdf', 'Path', path, 'dpi', 150, 'papersize', [l, w], 'Xsplit', 4,'Ysplit',3);

hold off

%% plot control action
figure
hold on
plot([0 20],[0 0],'k-.','LineWidth',1)
plot(t,uu,'color',[0 0.45 0.74])
plot([0 20],[-10 -10],'r--')
xlabel('t')
ylabel('$\tau$')
set(get(gca,'ylabel'),'rotation',0)
title('\textbf{Control action}')
legend('zero','control action','limit','Orientation','horizontal')
f= get(gcf);
f2p('con_a', 'Xlim', [0, 20], 'Ylim', [-10 2], 'Ytol', 0.05, 'Xtol', 0,...
        'extension', 'pdf', 'Path', path, 'dpi', 150, 'papersize', [l, w], 'Xsplit', 4,'Ysplit',3);

hold off