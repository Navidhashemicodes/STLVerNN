clear all
clc
close all

l0=[0.35; -0.35; 0.35];
u0=[0.4;  -0.3 ; 0.4];
n=3;
m=1;
timestep=0.2;
T=50;
controller_nn =  NN_Reader( 3,1,'neural_network_controller');


load('model_2_10.mat')
figure
for j=1:10000
    S2(:,1)=l0+rand(n,1).*(u0-l0);
    for i=1:T
        u=pred(controller_nn, S2(:,i));
        In=[S2(:,i);u];
        S2(:,i+1)=pred(Model_nn, In); 
    end
    
    plot3(S2(1,:),S2(2,:),S2(3,:))
    hold on
end
x=@(t,s) 0.1+0.2*t;
y=@(t,s) -(1.715+1.43*t)/14.3;
z=@(t,s) s;
s=fsurf(x,y,z,[0 1 -0.1 0.4]);%,'FaceAlpha',0.2);
% s.EdgeColor = 'none';
hold on

z=linspace(-0.1 ,  0.4  , 10);
y=linspace(-0.22, -0.12 , 10);
[Y, Z]=meshgrid(y,z);
X= 0.3*ones(size(Z));
s=surf(X,Y,Z);%,'FaceAlpha',0.2);
% s.EdgeColor = 'none';
hold on

x=linspace( 0.1, 0.3, 10);
z=linspace(-0.1, 0.4, 10);
[X, Z]=meshgrid(x,z);
Y= -0.12*ones(size(X));
s=surf(X,Y,Z);%,'FaceAlpha',0.2);
% s.EdgeColor = 'none';
hold on


x=@(t,s) 0.1+0.2*t;
y=@(t,s) -(1.715+1.43*t)/9.5;
z=@(t,s) s;
s=fsurf(x,y,z,[0 1 -0.1 0.4]);%,'FaceAlpha',0.2);
% s.EdgeColor = 'none';
hold on


z=linspace(-0.1 ,  0.4  , 10);
y=linspace(-0.33, -0.18 , 10);
[Y, Z]=meshgrid(y,z);
X= 0.1*ones(size(Z));
s=surf(X,Y,Z);%,'FaceAlpha',0.2);
% s.EdgeColor = 'none';
hold on

x=linspace( 0.1, 0.3, 10);
z=linspace(-0.1, 0.4, 10);
[X, Z]=meshgrid(x,z);
Y= -0.33*ones(size(X));
s=surf(X,Y,Z);%,'FaceAlpha',0.2);
% s.EdgeColor = 'none';
hold on



x=linspace( 0.02, 0.04 , 10);
y=linspace(-0.02, 0.02 , 10);
[X, Y]=meshgrid(x,y);
Z=-0.02*ones(size(X));
s=surf(X,Y,Z,'FaceAlpha',0.2);
s.EdgeColor = 'none';
hold on
Z=0.02*ones(size(X));
s=surf(X,Y,Z,'FaceAlpha',0.2);
s.EdgeColor = 'none';


x=linspace( 0.02, 0.04 , 10);
z=linspace(-0.02, 0.02 , 10);
[X, Z]=meshgrid(x,z);
Y=-0.02*ones(size(X));
s=surf(X,Y,Z,'FaceAlpha',0.2);
s.EdgeColor = 'none';
hold on
Y=0.02*ones(size(X));
s=surf(X,Y,Z,'FaceAlpha',0.2);
s.EdgeColor = 'none';

y=linspace(-0.02, 0.02 , 10);
z=linspace(-0.02, 0.02 , 10);
[Y, Z]=meshgrid(y,z);
X= 0.02*ones(size(X));
s=surf(X,Y,Z,'FaceAlpha',0.2);
s.EdgeColor = 'none';
hold on
X=0.04*ones(size(X));
s=surf(X,Y,Z,'FaceAlpha',0.2);
s.EdgeColor = 'none';

% figure
% for j=1:1000
%     S1(:,1)=l0+rand(n,1).*(u0-l0);
%     for i=1:T
%         
%         u=pred(controller_nn, S1(:,i));
%         [~,in_out] =  ode45(@(t,x)system_cont(t,x,u),[0 timestep],S1(:,i));
%         S1(:,i+1)=in_out(end,:)';
%         
%     end
%     
%     plot3(S1(1,:),S1(2,:),S1(3,:))
%     hold on
% end
% x=@(t,s) 0.1+0.2*t;
% y=@(t,s) -(1.715+1.43*t)/14.3;
% z=@(t,s) s;
% s=fsurf(x,y,z,[0 1 -0.1 0.4],'FaceAlpha',0.2);
% % s.EdgeColor = 'none';
% hold on
% 
% z=linspace(-0.1 ,  0.4  , 10);
% y=linspace(-0.22, -0.12 , 10);
% [Y, Z]=meshgrid(y,z);
% X= 0.3*ones(size(Z));
% s=surf(X,Y,Z,'FaceAlpha',0.2);
% % s.EdgeColor = 'none';
% hold on
% 
% x=linspace( 0.1, 0.3, 10);
% z=linspace(-0.1, 0.4, 10);
% [X, Z]=meshgrid(x,z);
% Y= -0.12*ones(size(X));
% s=surf(X,Y,Z,'FaceAlpha',0.2);
% % s.EdgeColor = 'none';
% hold on
% 
% 
% x=@(t,s) 0.1+0.2*t;
% y=@(t,s) -(1.715+1.43*t)/9.5;
% z=@(t,s) s;
% s=fsurf(x,y,z,[0 1 -0.1 0.4],'FaceAlpha',0.2);
% % s.EdgeColor = 'none';
% hold on
% 
% 
% z=linspace(-0.1 ,  0.4  , 10);
% y=linspace(-0.33, -0.18 , 10);
% [Y, Z]=meshgrid(y,z);
% X= 0.1*ones(size(Z));
% s=surf(X,Y,Z,'FaceAlpha',0.2);
% % s.EdgeColor = 'none';
% hold on
% 
% x=linspace( 0.1, 0.3, 10);
% z=linspace(-0.1, 0.4, 10);
% [X, Z]=meshgrid(x,z);
% Y= -0.33*ones(size(X));
% s=surf(X,Y,Z,'FaceAlpha',0.2);
% % s.EdgeColor = 'none';
% hold on
% 
% 
% 
% x=linspace( 0.02, 0.04 , 10);
% y=linspace(-0.02, 0.02 , 10);
% [X, Y]=meshgrid(x,y);
% Z=-0.02*ones(size(X));
% s=surf(X,Y,Z,'FaceAlpha',0.2);
% s.EdgeColor = 'none';
% hold on
% Z=0.02*ones(size(X));
% s=surf(X,Y,Z,'FaceAlpha',0.2);
% s.EdgeColor = 'none';
% 
% 
% x=linspace( 0.02, 0.04 , 10);
% z=linspace(-0.02, 0.02 , 10);
% [X, Z]=meshgrid(x,z);
% Y=-0.02*ones(size(X));
% s=surf(X,Y,Z,'FaceAlpha',0.2);
% s.EdgeColor = 'none';
% hold on
% Y=0.02*ones(size(X));
% s=surf(X,Y,Z,'FaceAlpha',0.2);
% s.EdgeColor = 'none';
% 
% y=linspace(-0.02, 0.02 , 10);
% z=linspace(-0.02, 0.02 , 10);
% [Y, Z]=meshgrid(y,z);
% X= 0.02*ones(size(X));
% s=surf(X,Y,Z,'FaceAlpha',0.2);
% s.EdgeColor = 'none';
% hold on
% X=0.04*ones(size(X));
% s=surf(X,Y,Z,'FaceAlpha',0.2);
% s.EdgeColor = 'none';

axis equal


function y = pred(net, x)
   
    len=length(net.weights)-1;
    for i=1:len
        x=poslin(net.weights{i}*x+net.biases{i});
    end
    y=net.weights{end}*x+net.biases{end};
    
end