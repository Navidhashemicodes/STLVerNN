function dxdt = system_cont(t,x,u)
    
    dxdt = [  -x(1)^3 + x(2); x(2)^3 + x(3); u];

end