clear all
close all
clc

N=1000;
M=100;
K=100;

% Sc1=1.001;
% Sc2=1;
% A=2*rand(2,2)-1;
% A=A/(Sc1*max(abs(eig(A))));
% B=Sc2*(2*rand(2,2)-1);
load('final_example.mat')
B=100*B;
figure
for i=1:N
%     thet=2*pi*i/N;
%     S(:,1)=K*[cos(thet) ; sin(thet)];
    thet=2*pi/3;
    S(:,1)=K*[cos(thet) ; sin(thet)]+10*rand(2,1);
    for j=1:M
        S(:,j+1)=A*S(:,j)+B*sign(S(:,j));
    end
    plot(S(1,:),S(2,:))
    hold on
    x(:,i)=S(:,1);
end
plot(x(1,:),x(2,:),'*')
axis equal