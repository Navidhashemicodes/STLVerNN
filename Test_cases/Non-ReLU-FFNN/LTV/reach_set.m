clear
clc
close all

load('LTV_disc.mat')


xxs =-1+1*rand(1,30);
yys = -1+1*rand(1,30);
num_steps=50;
t=0:T:num_steps*T;


% figure;
% hold on;
% 
% for xx = 1:length(xxs)
%     for yy = 1:length(yys)
%         s1(:,1) = [xxs(xx) yys(yy)];
%         for i=1:num_steps
%             u = sin(s1(1,i))-cos(s1(2,i));
%             [time,in_out] =  ode45(@(t,x)eq_fun(t+(i-1)*T,x,u),[0 T],s1(:,i)');
%             s1(:,i+1) = in_out(end,:)';
%         end
%         plot(t,s1(1,:),'r')
%         hold on
%         plot(t,s1(2,:),'b');
%         hold on; 
%     end
% end
% 
% 
% figure
% hold on
% for xx = 1:length(xxs)
%     for yy = 1:length(yys)
%         s(:,1) = [xxs(xx) yys(yy)];
%         for i=0:num_steps-1
%             u = sin(s(1,i+1))-cos(s(2,i+1));
%             j = mod(i,2*pi/T);
%             Phi = A_d{j+1};
%             B = B_d{j+1};
%             s(:,i+2)=Phi*s(:,i+1) + B*u;
%         end
%         plot(t,s(1,:),'r')
%         hold on
%         plot(t,s(2,:),'b');
%         hold on; 
%     end
% end


figure
hold on
load('neuro_controller.mat')
for xx = 1:length(xxs)
    for yy = 1:length(yys)
        s2(:,1) = [xxs(xx) yys(yy)];
        for i=0:num_steps-1
            u = pred_tanh(net,s2(:,i+1));
            j = mod(i,2*pi/T);
            Phi = A_d{j+1};
            B = B_d{j+1};
            s2(:,i+2)=Phi*s2(:,i+1) + B*u;
        end
%         plot(t,s2(1,:),'r')
        plot(s2(1,:),'r')
        hold on
%         plot(t,s2(2,:),'b');
        plot(s2(2,:),'b');
        hold on;
%         plot(s2(1,:),s2(2,:));
%         hold on;
        
    end
end




function y = pred_tanh(net, x)
   
    len=length(net.weights)-1;
    for i=1:len
        x=tanh(net.weights{i}*x+net.biases{i});
    end
    y=net.weights{end}*x+net.biases{end};
    
end