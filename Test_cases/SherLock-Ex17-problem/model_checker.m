clear all
clc
close all

l0=[0.35; -0.35;  0.35];
u0=[0.4 ; -0.3 ;  0.4];
n=3;
m=1;
S1(:,1)=l0+rand(3,1).*(u0-l0);
S2(:,1)=S1(:,1);
timestep=0.2;
T=50;
controller_nn =  NN_Reader( 3,1,'neural_network_controller');


model_type='normalized';
% model_type='regular';

switch model_type
    
    case 'regular'
        
%         load('Model.mat')
%         Model_nn.weights=net.weights;
%         Model_nn.biases=net.biases;
        load('networks.mat')
         
        for i=1:T
            
            u=pred(controller_nn, S1(:,i));
            [~,in_out] =  ode45(@(t,x)system_cont(t,x,u),[0 timestep],S1(:,i));
            S1(:,i+1)=in_out(end,:)';
            
            u=pred(controller_nn, S2(:,i));
            In=[S2(:,i);u];
            S2(:,i+1)=pred(Model_nn, In);
            
        end

        
        
        
    case 'normalized'
        
        load('Model_n.mat')
        Model_nn.weights=net.weights;
        Model_nn.biases=net.biases;
        
        load('maxmin.mat')
        a=-1; b=1;
        for i=1:T
            u=pred(controller_nn, S1(:,i));
            [~,in_out] =  ode45(@(t,x)system_cont(t,x,u),[0 timestep],S1(:,i));
            S1(:,i+1)=in_out(end,:)';
            
            u=pred(controller_nn, S2(:,i));
            In=[S2(:,i);u];
            In_n= (b-a) * diag(1./ (maxin-minin) ) * ( In - minin )  + a ;
            Out_n=pred(Model_nn, In_n);
            S2(:,i+1)= diag(maxout-minout)*((Out_n-a)/(b-a)) + minout;
            
            
        end
        Model_nn.biases{1}= -(b-a)*Model_nn.weights{1}*diag(1./ (maxin-minin) )*minin + a*Model_nn.weights{1}*ones(n+m,1) +Model_nn.biases{1};
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
save('model_2_10.mat', 'Model_nn')

for i=1:3
figure
plot(S1(i,:), 'r-')
hold on
plot(S2(i,:), 'blue-')
end
figure
plot(S2(1,:),S2(2,:))
figure
plot(S1(1,:),S1(2,:))
figure
plot(S2(1,:),S2(3,:))
figure
plot(S1(1,:),S1(3,:))
figure
plot(S2(2,:),S2(3,:))
figure
plot(S1(2,:),S1(3,:))
    
function y = pred(net, x)
   
    len=length(net.weights)-1;
    for i=1:len
        x=poslin(net.weights{i}*x+net.biases{i});
    end
    y=net.weights{end}*x+net.biases{end};
    
end