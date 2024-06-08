clear all
clc
close all

l0=[1;1];
u0=[2;2];

T=20;

load('Toy_example_2.mat')
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
%     plot(S2(1,5),S2(2,5), '*')
    hold on
end

x0=-0.8;y0=-0.5;x1=-0.4;y1=-0.7;
t=x0:0.001:x1;
plot(t, ((y1-y0)/(x1-x0))*(t-x0)+y0);

hold on

x0=-0.4;y0=-0.7;x1= 0.3;y1=-0.4;
t=x0:0.001:x1;
plot(t, ((y1-y0)/(x1-x0))*(t-x0)+y0);

hold on


x0=-0.1;y0=-0.2;x1= 0.3;y1=-0.4;
t=x0:0.001:x1;
plot(t, ((y1-y0)/(x1-x0))*(t-x0)+y0);

hold on


x0=-0.8;y0=-0.5;x1= -0.1;y1=-0.2;
t=x0:0.001:x1;
plot(t, ((y1-y0)/(x1-x0))*(t-x0)+y0);

hold on


 
plot(-1.1*ones(100,1), linspace(-0.3,0,100) )
hold on
plot(-0.7*ones(100,1), linspace(-0.3,0,100) )
hold on
plot(linspace(-1.1,-0.7,100), -0.3*ones(100,1)  )
hold on
plot(linspace(-1.1,-0.7,100),  0*ones(100,1)  )

axis equal
xlim([-1.5 2.5])
ylim([-1.5 2.5])

grid on

ax = gca;

ax.LineWidth = 2;
ax.FontSize = 18;
ax.XTick = [-1.5 -1 -0.5  0    0.5  1   1.5  2  2.5];
ax.YTick = [-1.5 -1 -0.5  0    0.5  1   1.5  2  2.5];


function y = pred(net, x)
   
    len=length(net.weights)-1;
    for i=1:len
        x=tanh(net.weights{i}*x+net.biases{i});
    end
    y=net.weights{end}*x+net.biases{end};
    
end