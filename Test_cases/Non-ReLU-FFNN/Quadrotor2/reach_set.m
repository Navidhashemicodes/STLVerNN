clear all
clc
close all



l0=[0.0638; 0.0638;  -0.0231; 0;  0;  0];
u0=[0.1063; 0.1063 ;  0.0231; 0;  0;  0];
n=6;
m=3;

timestep=0.05;
T=20;
FF=0.1*0.5*0.85;
load('networks.mat')

figure
for j=1:10000
    S2(:,1)=l0+rand(n,1).*(u0-l0);
    for i=1:T
        u=pred(nn(i), S2(:,i));
        In=[S2(:,i);u];
        S2(:,i+1)=pred(Model_nn, In);
    end
    
    plot3(S2(1,:),S2(2,:),S2(3,:))
    hold on
end

plotcube(FF*[2;8;0], FF*2.5 , 'blue' , 0.2);

hold on

plotcube(FF*[8;2;0], FF*2.5 , 'blue' , 0.2);

% hold on
% 
% plotcube(FF*[8;8;-3], FF*2.5 , 'green' , 0.2);

hold on

plotcube(FF*[5;5;0], FF*2.5 , 'red' , 0.2);

hold on

plotcube(FF*[2;2;0], FF*1 , 'black' , 0.2);


xlim([0,0.4])
ylim([0,0.4])
zlim([-0.2,0.2])


% figure
% for j=1:1000
%     S1(:,1)=l0+rand(n,1).*(u0-l0);
%     for i=1:T
%         
%         u=pred(nn(i), S1(:,i));
%         [~,in_out] =  ode45(@(t,x)system_cont(t,x,u),[0 timestep],S1(:,i));
%         S1(:,i+1)=in_out(end,:)';
%         
%     end
%     
%     plot3(S1(1,:),S1(2,:),S1(3,:))
%     hold on
% end
% 
% plotcube(FF*[2;8;0], FF*2.5 , 'blue' , 0.2);
% 
% hold on
% 
% plotcube(FF*[8;2;0], FF*2.5 , 'blue' , 0.2);
% 
% % hold on
% % 
% % plotcube(FF*[8;8;-3], FF*2.5 , 'green' , 0.2);
% 
% hold on
% 
% plotcube(FF*[5;5;0], FF*2.5 , 'red' , 0.2);
% 
% hold on
% 
% plotcube(FF*[2;2;0], FF*1 , 'black' , 0.2);
% 
% axis equal


function y = pred(net, x)
   
    len=length(net.weights)-1;
    for i=1:len
        x=tanh(net.weights{i}*x+net.biases{i});
    end
    y=net.weights{end}*x+net.biases{end};
    
end