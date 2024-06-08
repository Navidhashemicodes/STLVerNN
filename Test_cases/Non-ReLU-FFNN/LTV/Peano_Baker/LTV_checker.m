clear
clc
close all
x0 = 10*rand(2,1);
horiz=10*pi;

% [tout,xout] = ode45(@(t,x)eq_fun(t,x,0),[0, horiz],x0');

load('A_d.mat')   %%% T is also here
% n=4;
% T=2*pi/20;
% N=150;
% parfor i=1:2*pi/T
%     A_d{i} = peano_baker(i-1, T, N, n);
% end


num_steps=horiz/T;

s1(:,1)=x0;
for i=1:num_steps
%     u=sin(s1(1,i))-cos(s1(2,i));
    [time,in_out] =  ode45(@(t,x)eq_fun(t+(i-1)*T,x,0),[0 T],s1(:,i)');
    s1(:,i+1) = in_out(end,:)';
end



s(:,1)=x0;
for i=0:num_steps-1
    j=mod(i,2*pi/T);
    Phi=A_d{j+1};
%     u=sin(s(1,i+1))-cos(s(2,i+1));
    s(:,i+2)=Phi*s(:,i+1) + T*[1;0]*0;
end

t=0:T:horiz;

figure
% plot(tout,xout(:,1),'black')
% hold on
plot(t,s1(1,:),'r')
hold on
plot(t,s(1,:), 'b')


figure
% plot(tout,xout(:,2),'black')
% hold on
plot(t,s1(2,:),'r')
hold on
plot(t,s(2,:), 'b')
