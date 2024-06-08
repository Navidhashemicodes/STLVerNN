clear all
close all
clc

N=1000;
M=100;

load('final_example_BigB.mat')
load('network.mat')

figure
for i=1:N
    S(:,1)=[-45;90]+(2*rand(2,1)-1).*[5;5];
    for j=1:M
        S(:,j+1)=A*S(:,j)+B*NN(nn_controller,S(:,j));
    end
    plot(S(1,:),S(2,:))
    hold on
    x(:,i)=S(:,1);
end
% plot(x(1,:),x(2,:),'*')
hold on


%%%%  P1
b1=-250;a1=35;b2=180;a2=-45;
x=0.3*linspace(b1-a1,b2-a1,50);
y=x/3+a1;
plot(x,y, 'black');
hold on
x=0.3*linspace(b1-a2,b2-a2,50);
y=x/3+a2;
plot(x,y, 'black');
hold on

x=0.3*linspace(b1-a1,b1-a2,50);
y=-3*x+b1;
plot(x,y, 'black');
hold on
x=0.3*linspace(b2-a1,b2-a2,50);
y=-3*x+b2;
plot(x,y, 'black');
hold on



%%%%  P2
b1=-130;a1=16;b2=80;a2=-22;
x=0.3*linspace(b1-a1,b2-a1,50);
y=x/3+a1;
plot(x,y, 'black');
hold on
x=0.3*linspace(b1-a2,b2-a2,50);
y=x/3+a2;
plot(x,y, 'black');
hold on

x=0.3*linspace(b1-a1,b1-a2,50);
y=-3*x+b1;
plot(x,y, 'black');
hold on
x=0.3*linspace(b2-a1,b2-a2,50);
y=-3*x+b2;
plot(x,y, 'black');
hold on



%%%%  P3
b1=-55;a1=7;b2=40;a2=-10;
x=0.3*linspace(b1-a1,b2-a1,50);
y=x/3+a1;
plot(x,y, 'black');
hold on
x=0.3*linspace(b1-a2,b2-a2,50);
y=x/3+a2;
plot(x,y, 'black');
hold on

x=0.3*linspace(b1-a1,b1-a2,50);
y=-3*x+b1;
plot(x,y, 'black');
hold on
x=0.3*linspace(b2-a1,b2-a2,50);
y=-3*x+b2;
plot(x,y, 'black');
hold on



%%%%  P4
b1=-25;a1=3;b2=15;a2=-4;
x=0.3*linspace(b1-a1,b2-a1,50);
y=x/3+a1;
plot(x,y, 'black');
hold on
x=0.3*linspace(b1-a2,b2-a2,50);
y=x/3+a2;
plot(x,y, 'black');
hold on

x=0.3*linspace(b1-a1,b1-a2,50);
y=-3*x+b1;
plot(x,y, 'black');
hold on
x=0.3*linspace(b2-a1,b2-a2,50);
y=-3*x+b2;
plot(x,y, 'black');
hold on




%%%%  P5
b1=-9;a1=1.1;b2=5;a2=-1.2;
x=0.3*linspace(b1-a1,b2-a1,50);
y=x/3+a1;
plot(x,y, 'black');
hold on
x=0.3*linspace(b1-a2,b2-a2,50);
y=x/3+a2;
plot(x,y, 'black');
hold on

x=0.3*linspace(b1-a1,b1-a2,50);
y=-3*x+b1;
plot(x,y, 'black');
hold on
x=0.3*linspace(b2-a1,b2-a2,50);
y=-3*x+b2;
plot(x,y, 'black');
hold on



%%%%  P6
b1=-1.2;a1=0.7;b2=3.5;a2=-0.6;
x=0.3*linspace(b1-a1,b2-a1,50);
y=x/3+a1;
plot(x,y, 'black');
hold on
x=0.3*linspace(b1-a2,b2-a2,50);
y=x/3+a2;
plot(x,y, 'black');
hold on

x=0.3*linspace(b1-a1,b1-a2,50);
y=-3*x+b1;
plot(x,y, 'black');
hold on
x=0.3*linspace(b2-a1,b2-a2,50);
y=-3*x+b2;
plot(x,y, 'black');
hold on




axis equal
ax= gca;
xlim([-150,125])
ylim([-100,100])
ax.YTick =[-100 -75  -50   -25   0   25  50   75   100];
ax.XTick =[-150, -125, -100, -75 , -50, -25, 0 , 25, 50, 75, 100, 125 ];
grid on

ax.LineWidth = 3;
ax.FontSize = 18;
