clear all
clc
close all

l0=[90;32;0;10;30;0];
u0=[110;32.2;0;11;30.2;0];
S1(:,1)=l0+rand(6,1).*(u0-l0);
S2(:,1)=S1(:,1);
timestep=0.1;
T=50;
V_set=30;
t_gap=1.4;
E=[ 0  0  0  0  1  0;...
    1  0  0 -1  0  0;...
    0  1  0  0 -1  0];


load('controller_3_20.mat')
controller_nn.weights=network.W;
controller_nn.biases=network.b;


% model_type='normalized';
model_type='regular';

switch model_type
    
    case 'regular'
        
%         load('Model.mat')
%         Model_nn.weights=net.weights;
%         Model_nn.biases=net.biases;
        load('networks.mat')
         
        for i=1:T
            init_a=[V_set; t_gap; E*S1(:,i)];
            a_ego=pred(controller_nn, init_a);
            [~,in_out] =  ode45(@(t,x)dynamicsACC(t,x,a_ego),[0 timestep],S1(:,i));
            S1(:,i+1)=in_out(end,:)';
            
            init_a=[V_set; t_gap; E*S2(:,i)];
            a_ego=pred(controller_nn, init_a);
            In=[S2(:,i);a_ego];
            S2(:,i+1)=pred(Model_nn, In);
        end

        
        
        
    case 'normalized'
        
        load('Model_n.mat')
        Model_nn.weights=net.weights;
        Model_nn.biases=net.biases;
        
        load('maxmin.mat')
        a=-1; b=1;
        for i=1:T
            init_a=[V_set; t_gap; E*S1(:,i)];
            a_ego=pred(controller_nn, init_a);
            [~,in_out] =  ode45(@(t,x)dynamicsACC(t,x,a_ego),[0 timestep],S1(:,i));
            S1(:,i+1)=in_out(end,:)';
            
            init_a=[V_set; t_gap; E*S2(:,i)];
            a_ego=pred(controller_nn, init_a);
            In=[S2(:,i);a_ego];
            In_n= (b-a) * diag(1./ (maxin-minin) ) * ( In - minin )  + a ;
            Out_n=pred(Model_nn, In_n);
            S2(:,i+1)= diag(maxout-minout)*((Out_n-a)/(b-a)) + minout;
            
            
        end
        Model_nn.biases{1}= -(b-a)*Model_nn.weights{1}*diag(1./ (maxin-minin) )*minin + a*Model_nn.weights{1}*ones(7,1) +Model_nn.biases{1};
        Model_nn.weights{1}= (b-a)*Model_nn.weights{1}*diag(1./ (maxin-minin) );
        Model_nn.weights{end}= diag(maxout-minout)*Model_nn.weights{end}/(b-a);
        Model_nn.biases{end}=  diag(maxout-minout)*(Model_nn.biases{end}-a)/(b-a) + minout;
        
        
end

for i=1:length(Model_nn.biases)-1
    funcs=cell(size(Model_nn.biases{i},1),1);
    funcs(:)={'poslin'};
    Model_nn.layers{i}=funcs;
end
for i=1:length(controller_nn.biases)-1
    funcs=cell(size(controller_nn.biases{i},1),1);
    funcs(:)={'poslin'};
    controller_nn.layers{i}=funcs;
end
save('networks.mat', 'Model_nn', 'controller_nn')

for i=1:6
figure
plot(S1(i,:), 'r-')
hold on
plot(S2(i,:), 'blue-')
end
    
function y = pred(net, x)
   
    len=length(net.weights)-1;
    for i=1:len
        x=poslin(net.weights{i}*x+net.biases{i});
    end
    y=net.weights{end}*x+net.biases{end};
    
end