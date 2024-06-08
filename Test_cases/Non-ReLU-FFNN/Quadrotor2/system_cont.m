function dxdt = system_cont(t,x,u)
    
    dxdt = [  x(4) ;  x(5) ;  x(6) ; 9.81*tan(u(1)) ; -9.81*tan(u(2)) ; u(3)-9.81];

end