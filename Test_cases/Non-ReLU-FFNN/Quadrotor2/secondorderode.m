syms x(t) u  x10  x20
Dx = diff(x);

ode = diff(Dx,t) == Dx-(2+sin(t))*x;
cond1 = x(0) == x10;
cond2 = Dx(0) == x20+u;


conds = [cond1 cond2];
sol=dsolve(ode,'Implicit',true)
% xSol(t) = dsolve(ode,conds);
% xSol = simplify(xSol)