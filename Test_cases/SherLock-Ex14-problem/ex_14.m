function dxdt = ex_14(t,x,u)
 
    dxdt =[ -x(1)*(0.1 + (x(1) + x(2))^2) ; (u + x(1)) * (0.1 + (x(1) + x(2))^2) ];

end