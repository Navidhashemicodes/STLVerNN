function dxdt = eq_fun(t,x,u)
    dxdt = [x(2)+u  ;  -x(2)-(2+sin(t))*x(1)];
end