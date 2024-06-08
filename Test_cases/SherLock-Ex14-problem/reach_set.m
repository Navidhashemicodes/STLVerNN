clear all
clc
% close all

l0=[0.8;0.4];
u0=[0.9;0.5];

timestep=0.1;
T=100;
controller_nn =  NN_Reader( 2,1,'neural_network_controller');


load('model_3_10.mat')

for j=1:5000
    S2(:,1)=l0+rand(2,1).*(u0-l0);
    for i=1:T
        u=pred(controller_nn, S2(:,i));
        In=[S2(:,i);u];
        S2(:,i+1)=pred(Model_nn, In);
    end
    
    plot(S2(1,:),S2(2,:));
    hold on
end

t=0:0.01:1;
plot(0.5+0.15*t, 0.3*t )
hold on
plot(0.5*ones(100,1), linspace(0,0.3,100) )
hold on
plot(linspace(0.5,0.65,100), 0.3*ones(100,1)  )


hold on
plot(0.7+0.1*t , 0.1*t  )
hold on
plot(linspace(0.7,0.8,100) , zeros(100,1) )
hold on
plot(0.8*ones(100,1), linspace(0,0.1,100) )


hold on
plot(linspace(0.1,0.15,100), -0.02*ones(100,1)  )
hold on
plot(linspace(0.1,0.15,100),  0.02*ones(100,1)  )
hold on
plot(0.1*ones(100,1) , linspace(-0.02,0.02,100))
hold on
plot(0.15*ones(100,1) , linspace(-0.02,0.02,100))

axis equal
xlim([-0.1,1  ])
ylim([-0.2,0.6])

axis equal
xlim([0 1])
ylim([-0.4 0.6])

grid on

ax = gca;

ax.LineWidth = 4;
ax.FontSize = 24;
ax.XTick = [ 0    0.2  0.4   0.6  0.8  1     ];
ax.YTick = [-0.4 -0.2  0     0.2  0.4  0.6   ];


% figure
% for j=1:1000
%     S1(:,1)=l0+rand(2,1).*(u0-l0);
%     for i=1:T
%         
%         u=pred(controller_nn, S1(:,i));
%         [~,in_out] =  ode45(@(t,x)ex_14(t,x,u),[0 timestep],S1(:,i));
%         S1(:,i+1)=in_out(end,:)';
%         plot(S1(1,i),S1(2,i), char(col(mod(i,5)+1)));
%         hold on
%     end
%      
% %     plot(S1(1,:),S1(2,:),char(col(mod(j,5)+1)))
% %     hold on
% end
% t=0:0.01:1;
% plot(0.5+0.15*t, 0.3*t )
% hold on
% plot(0.5*ones(100,1), linspace(0,0.3,100) )
% hold on
% plot(linspace(0.5,0.65,100), 0.3*ones(100,1)  )
% 
% 
% hold on
% plot(0.7+0.1*t , 0.1*t  )
% hold on
% plot(linspace(0.7,0.8,100) , zeros(100,1) )
% hold on
% plot(0.8*ones(100,1), linspace(0,0.1,100) )
% 
% 
% 
% hold on
% plot(linspace(0.1,0.15,100), -0.02*ones(100,1)  )
% hold on
% plot(linspace(0.1,0.15,100),  0.02*ones(100,1)  )
% hold on
% plot(0.1*ones(100,1) , linspace(-0.02,0.02,100))
% hold on
% plot(0.15*ones(100,1) , linspace(-0.02,0.02,100))
% 
% axis equal
% xlim([-0.1,1  ])
% ylim([-0.2,0.6])



function y = pred(net, x)
   
    len=length(net.weights)-1;
    for i=1:len
        x=poslin(net.weights{i}*x+net.biases{i});
    end
    y=net.weights{end}*x+net.biases{end};
    
end