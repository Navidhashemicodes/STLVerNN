clear all
clc
close all

l0=[1;1];
u0=[2;2];

T=30;

load('Best_example.mat')
controller_nn.weights=NN_actor.weights;
Model_nn.weights=NN_model.weights;
controller_nn.biases=NN_actor.Biases;
Model_nn.biases=NN_model.Biases;

figure
for j=1:10000
    S2(:,1)=l0+rand(2,1).*(u0-l0);
    for i=1:T
        u=pred(controller_nn, S2(:,i));
        In=[S2(:,i);u];
        S2(:,i+1)=pred(Model_nn, In); 
    end
    
    plot(S2(1,:),S2(2,:))
    hold on
end


%%%P1
plot(-2*ones(100,1), linspace(-0.8,-0.2,100) )
hold on
plot(-4*ones(100,1), linspace(-0.8,-0.2,100) )
hold on
plot(linspace(-4,-2,100), -0.8*ones(100,1)  )
hold on
plot(linspace(-4,-2,100),  -0.2*ones(100,1)  )
hold on

%%%P2 
plot(-7*ones(100,1), linspace(-0.8,-0.2,100) )
hold on
plot(-10*ones(100,1), linspace(-0.8,-0.2,100) )
hold on
plot(linspace(-10,-7,100), -0.8*ones(100,1)  )
hold on
plot(linspace(-10,-7,100),  -0.2*ones(100,1)  )
hold on


%%%P3
plot(-14*ones(100,1), linspace(-0.6,-0.2,100) )
hold on
plot(-24*ones(100,1), linspace(-0.6,-0.2,100) )
hold on
plot(linspace(-24,-14,100), -0.6*ones(100,1)  )
hold on
plot(linspace(-24,-14,100),  -0.2*ones(100,1)  )
hold on


% plot(-32*ones(100,1), linspace(-2,1,100) )
% hold on
% plot(-42*ones(100,1), linspace(-2,1,100) )
% hold on
% plot(linspace(-42,-32,100), -2*ones(100,1)  )
% hold on
% plot(linspace(-42,-32,100),  1*ones(100,1)  )
% hold on
% 
% 
%%%P5
plot(-60*ones(100,1), linspace(0.2,0.6,100) )
hold on
plot(-75*ones(100,1), linspace(0.2,0.6,100) )
hold on
plot(linspace(-75,-60,100), 0.2*ones(100,1)  )
hold on
plot(linspace(-75,-60,100),  0.6*ones(100,1)  )
hold on


% axis equal

ax= gca
ax.XScale='linear'

xlim([-78,3])
ylim([-2,4])
grid on
ax.LineWidth = 2;
ax.FontSize = 18;
% ax.XTick = [-1.5 -1 -0.5 0    0.5 1   1.5];
% ax.YTick = [-1 -0.5  0   0.5  1   1.5 2 ];

function y = pred(net, x)
   
    len=length(net.weights)-1;
    for i=1:len
        x=poslin(net.weights{i}*x+net.biases{i});
    end
    y=net.weights{end}*x+net.biases{end};
    
end