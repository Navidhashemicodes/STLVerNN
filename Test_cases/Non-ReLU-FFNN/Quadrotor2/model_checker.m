clear all
clc
close all


l0=[0.0638; 0.0638;  -0.0231; 0;  0;  0];
u0=[0.1063; 0.1063 ;  0.0231; 0;  0;  0];
n=6;
m=3;
S1(:,1)=l0+rand(6,1).*(u0-l0);
S2(:,1)=S1(:,1);
timestep=0.05;
T=55;



load('nn_controller.mat')

for i=1:55
    NN_actors=NN_actor;
    weights1 = NN_actors.weights{1};
    biases1 = NN_actors.biases{1};
    biases1 = biases1+(i-1)*weights1(:,7);
    weights1 = weights1(:,1:6);
    NN_actors.weights{1} = weights1;
    NN_actors.biases{1} = biases1;
    nn(i).weights = NN_actors.weights;
    nn(i).biases = NN_actors.biases;
    nn(i).layers = NN_actors.layers;
end




model_type='normalized';
% model_type='regular';

switch model_type
    
    case 'regular'
        
        load('Model.mat')
        Model_nn.weights=net.weights;
        Model_nn.biases=net.biases;
%         load('networks.mat')
         
        for i=1:T
            
            u=pred(nn(i), S1(:,i));
            [~,in_out] =  ode45(@(t,x)system_cont(t,x,u),[0 timestep],S1(:,i));
            S1(:,i+1)=in_out(end,:)';
            
            u=pred(nn(i), S2(:,i));
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
            u=pred(nn(i), S1(:,i));
            [~,in_out] =  ode45(@(t,x)system_cont(t,x,u),[0 timestep],S1(:,i));
            S1(:,i+1)=in_out(end,:)';
            
            u=pred(nn(i), S2(:,i));
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
    funcs(:)={'tanh'};
    Model_nn.layers{i}=funcs;
end


save('networks.mat', 'Model_nn', 'nn')

for i=1:6
figure
plot(S1(i,:), 'r-')
hold on
plot(S2(i,:), 'blue-')
end
figure
plot3(S2(1,:),S2(2,:),S2(3,:), 'black')
hold on
plot3(S1(1,:),S1(2,:),S1(3,:), 'blue')

    
function y = pred(net, x)
   
    len=length(net.weights)-1;
    for i=1:len
        x=tanh(net.weights{i}*x+net.biases{i});
    end
    y=net.weights{end}*x+net.biases{end};
    
end