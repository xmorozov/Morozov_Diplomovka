%arm 9
%darm 5 'Ylim', [-300, 200]
%dpend 'Ylim', [-1200, 800] 5
f2p('control','Ylim', [-10, 10], 'Ytol', 0.05, 'Xtol', 0,...
'extension', 'pdf', 'dpi', 150, 'Path', path, 'papersize', [l, w], 'Xsplit', 10,'Ysplit',4);
