function dydt = odefcn(t,y,A,B,x)
dydt = zeros(2,1);
dydt(1) = y(2)+x;
dydt(2) = (A/B)*t.*y(1);
end