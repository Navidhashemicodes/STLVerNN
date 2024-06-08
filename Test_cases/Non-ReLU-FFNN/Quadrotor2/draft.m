% clear all
% clc
% close all
% 
% load('results.mat')
% clearvars -except param
% n=6;
% M=10;
% P=3;
% Param=param{end};
% 
% index=0;
% Weights{1}=reshape(Param(1,index+1:index+M*(n+1)), [n+1, M])';  
% index=index+M*(n+1);
% Biases{1}= Param(1,index+1:index+M)';
% index=index+M;
% Weights{2}=reshape(Param(1,index+1:index+M*P), [M, P])';
% Weights{2}=0.1*eye(3)*Weights{2};
% index=index+M*P;
% Biases{2} = Param(1,index+1:index+P)';
% Biases{2} = 0.1*eye(3)*Biases{2};
% Weights{3} = diag([0.1,0.1,2]);
% Biases{3} = [0; 0 ; 9.81];
% 
% Layers{1}=cell(10,1);
% Layers{1}(:)={'tanh'};
% 
% Layers{2}=cell(3,1);
% Layers{2}(:)={'tanh'};
% 
% 
% NN_actor.weights=Weights;
% NN_actor.biases=Biases;
% NN_actor.layers=Layers;
% 
% 
% save('nn_controller.mat', 'NN_actor')
% clear all
% clc
% close all


clear
clc

hh=64;

n1=4;
init_dim=3;

hh=hh-1;
rr=dec2base(hh,n1);
len=length(rr);
R=zeros(1,init_dim);
if len==1
    R(1)=0;
    R(2)=0;
    R(3)=str2num(rr(1));
elseif len==2
    R(1)=0;
    R(2)=str2num(rr(1));
    R(3)=str2num(rr(2));
elseif len==3
    R(1)=str2num(rr(1));
    R(2)=str2num(rr(2));
    R(3)=str2num(rr(3));
end
v=R+1