function A_d = peano_baker(k, T, N, n)
%%% N  is number of segmentation on T 
dt= T/N;

A_d = recursive_integrator(k, N+1, 0, T, dt);

for i=1:n
    A_d = A_d + recursive_integrator(k, N+1, i, T, dt);
end

end


function  A_delta = recursive_integrator(k,k_p_order,order,T, dt)

if order==0
    A_delta=eye(2);
    return;
end

tau=k*T + (0:dt:k_p_order*dt)  ;

A_delta=[0  1 ; -(2+sin(tau(1)))   -1]   *  recursive_integrator( k , (tau(1)-k*T)/dt, order-1, T, dt )    *   dt;
for i=2:k_p_order+1
    A_delta=A_delta+[0  1 ; -(2+sin(tau(i)))   -1]   *   recursive_integrator( k , (tau(i)-k*T)/dt, order-1, T, dt )    *   dt;
end


end