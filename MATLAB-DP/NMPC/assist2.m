%arm 9
%darm 5 'Ylim', [-300, 200]
%dpend 'Ylim', [-600, 100] 7
%pend [-50, 200] 5
%control [-10, 10] 4
w = 9;
l = 42;

f2p('arm','Ylim', [-30, 30], 'Ytol', 0.05, 'Xtol', 0,...
'extension', 'pdf', 'dpi', 150, 'Path', path, 'papersize', [l, w], 'Xsplit',10,'Ysplit',4);
%%
f2p('darm','Ylim', [-40, 120], 'Ytol', 0.05, 'Xtol', 0,...
'extension', 'pdf', 'dpi', 150, 'Path', path, 'papersize', [l, w], 'Xsplit',10,'Ysplit',8);
%%
f2p('pend','Ylim', [-50, 200], 'Ytol', 0.05, 'Xtol', 0,...
'extension', 'pdf', 'dpi', 150, 'Path', path, 'papersize', [l, w], 'Xsplit',10,'Ysplit',5);
%%
f2p('dpend','Ylim', [-140, 20], 'Ytol', 0.05, 'Xtol', 0,...
'extension', 'pdf', 'dpi', 150, 'Path', path, 'papersize', [l, w], 'Xsplit',10,'Ysplit',8);
%%
f2p('control','Ylim', [-10, 10], 'Ytol', 0.05, 'Xtol', 0,...
'extension', 'pdf', 'dpi', 150, 'Path', path, 'papersize', [l, w], 'Xsplit',10,'Ysplit',4);

%%
hold on
plot([0 15],[10 10],'k--')
plot([0 15],[-10 -10],'k--')