function dxdt = AxBsignx(t,A,B,x)
 
    dxdt =A*[x(1);x(2)]+B*sign([x(1);x(2)]);

end