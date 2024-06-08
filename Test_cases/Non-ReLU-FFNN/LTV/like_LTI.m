function [A_d, B_d] = like_LTI(k, N, T)
%%% N  is number of segmentation on T 
dt= T/N;

A_d = expm(integrator1(k,N+1, T , dt));

B_d = integrator2(k, N+1 ,T, dt)*[1;0];

end


function  A_delta = integrator1(k, k_p ,T, dt)

tau=k*T + (0:dt:k_p*dt)  ;

A_delta=[0  1 ; -(2+sin(tau(1)))   -1]    *   dt;
for i=2:k_p+1
    A_delta=A_delta+[0  1 ; -(2+sin(tau(i)))   -1]  *   dt;
end


end


function  B_delta = integrator2(k, k_p ,T, dt)


B_delta=expm(integrator1(k, 0, T, dt))    *   dt;
for i=1:k_p
    B_delta=B_delta + expm(integrator1(k, i, T, dt))    *   dt;
end


end