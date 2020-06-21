%arm 9 'Ylim', [-140, 40]
%darm 5 'Ylim', [-300, 200]
%dpend 'Ylim', [-1200, 800] 5
w = 9;
l = 42;

f2p('arm', 'Xlim', [0, tf],'Ylim', [-100, 40], 'Ytol', 0.05, 'Xtol', 0,...
'extension', 'pdf', 'dpi', 150, 'Path', path, 'papersize', [l, w], 'Xsplit',10,'Ysplit',7);
%%
f2p('darm', 'Xlim', [0, tf],'Ylim', [-40, 120], 'Ytol', 0.05, 'Xtol', 0,...
'extension', 'pdf', 'dpi', 150, 'Path', path, 'papersize', [l, w], 'Xsplit',10,'Ysplit',8);
%%
f2p('pend', 'Xlim', [0, tf],'Ylim', [-50, 300], 'Ytol', 0.05, 'Xtol', 0,...
'extension', 'pdf', 'dpi', 150, 'Path', path, 'papersize', [l, w], 'Xsplit',10,'Ysplit',7);
%%
f2p('dpend', 'Xlim', [0, tf],'Ylim', [-140, 20], 'Ytol', 0.05, 'Xtol', 0,...
'extension', 'pdf', 'dpi', 150, 'Path', path, 'papersize', [l, w], 'Xsplit',10,'Ysplit',8);
%%
f2p('control', 'Xlim', [0, tf],'Ylim', [-10, 4], 'Ytol', 0.05, 'Xtol', 0,...
'extension', 'pdf', 'dpi', 150, 'Path', path, 'papersize', [l, w], 'Xsplit',10,'Ysplit',4);

%%
w = 9;
l = 42;
x1 = 0.32
y1 = 157
x = 0.74
y = 248
hold on 
plot([0 4],[180 180],'k--')
plot([x x],[180 y],'r-','linewidth',3)
plot(x1,y1,'ro','markersize',6,'linewidth',4)
plot(x,y,'ro','markersize',6,'linewidth',4,'MarkerFaceColor','r')
f2p('swing3', 'Xlim', [0, tf],'Ylim', [-200, 500], 'Ytol', 0.05, 'Xtol', 0,...
'extension', 'pdf', 'dpi', 150, 'Path', path, 'papersize', [l, w], 'Xsplit',8,'Ysplit',4);

%%
hold on
plot([0 15],[10 10],'k--')
plot([0 15],[-10 -10],'k--')