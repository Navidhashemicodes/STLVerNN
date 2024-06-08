clear all
close all
clc

N=1;
M=1000;
K=10;

Sc1=10;
Sc2=0.01;
V=2*rand(2,2)-1;
E=-Sc1*diag( rand(2,1));
A=inv(V)*E*V;
B=Sc2*(2*rand(2,2)-1);
timestep=0.1;
figure
for i=1:N
    thet=2*pi*i/N;
    S(:,1)=K*[cos(thet) ; sin(thet)];
    for j=1:M
        [~,in_out] =  ode45(@(t,x)AxBsignx(t,A,B,x),[0 timestep],S(:,j)');
        S(:,j+1)=in_out(end,:)';
    end
    plot(S(1,:),S(2,:))
    hold on
    x(:,i)=S(:,1);
end
plot(x(1,:),x(2,:),'*')
axis equal