clear
clc
close all
x0 = 3*(2*rand(2,1)-1);
horiz=150*pi;



% T=2*pi/30;
% N=1000;
% 
% parfor i=1:2*pi/T
%     [A_d{i}, B_d{i}] = like_LTI(i-1, N+1, T);
% end
load('LTV_disc.mat')




num_steps=horiz/T;

s1(:,1)=x0;
for i=1:num_steps
    u = sin(s1(1,i))-cos(s1(2,i));
    [time,in_out] =  ode45(@(t,x)eq_fun(t+(i-1)*T,x,u),[0 T],s1(:,i)');
    s1(:,i+1) = in_out(end,:)';
end



s(:,1)=x0;
for i=0:num_steps-1
    u = sin(s(1,i+1))-cos(s(2,i+1));
    j = mod(i,2*pi/T);
    Phi = A_d{j+1};
    B = B_d{j+1};
    s(:,i+2)=Phi*s(:,i+1) + B*u;
end

load('neuro_controller.mat')

s2(:,1)=x0;
for i=0:num_steps-1
    u = pred_tanh(net, s2(:,i+1));
    j = mod(i,2*pi/T);
    Phi = A_d{j+1};
    B = B_d{j+1};
    s2(:,i+2)=Phi*s2(:,i+1) + B*u;
end

t=0:T:horiz;

figure
plot(t,s2(1,:),'black')
hold on
plot(t,s1(1,:),'r')
hold on
% plot(t,s(1,:), 'b')


figure
plot(t,s2(2,:),'black')
hold on
plot(t,s1(2,:),'r')
hold on
% plot(t,s(2,:), 'b')



function y = pred_tanh(net, x)
   
    len=length(net.weights)-1;
    for i=1:len
        x=tanh(net.weights{i}*x+net.biases{i});
    end
    y=net.weights{end}*x+net.biases{end};
    
end