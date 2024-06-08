clear all
clc
close all

for i=1:100
    model.A{i}=1;
    model.B{i}=0;
    model.E{i}=1;
    model.C{i}=1;
    model.D{i}=0;
end
model.type='LTV';

nncontroller.weights={1,1};
nncontroller.biases ={0,0};
func=cell(1,1);
func(:)={'poslin'};
nncontroller.layers{1}=func;

%%%%% be careful about the number of cells per parameter. It should be
%%%%% consistant with the STL specification you will introduce later,




str11= str2stl('[1,0]');
str= G_operation(str11, 0, 5);

p = '( []_[0,5]( p1 ) )';


predicate_transformation.C(1,:)= -1;   %%% x-0.1>0
predicate_transformation.d(1,1)= 10;


p_map(1).str= 'p1';
p_map(1).A  =   1;
p_map(1).b  =  10;




[net, border] = Trapezius_maker(model, nncontroller, predicate_transformation , str);
Net=Linear_Nonlinear(net, border);


for kk=1:1
s(:,1)=1;
for i=1:5
    a(:,i)=pred(nncontroller,s(:,i));
    s(:,i+1)=s(:,i)+1;
end

seqT=[0 , 1 , 2 , 3 , 4 , 5]';
seqS=[1 , 2 , 3 , 4 , 5 , 6]';
[rob,~] = dp_taliro(p,p_map,seqS,seqT);

G=predicate_transformation.C*s(:,1)+predicate_transformation.d;

for i=2:6
    G=min(G, predicate_transformation.C*s(:,i)+predicate_transformation.d);
end
clc
rho1(kk)=G;
rho2=NN(Net, s(:,1));
disp(rob)
disp(rho1(end))
disp(rho2)
end 
    
function y = pred(net, x)
   
    len=length(net.weights)-1;
    for i=1:len
        x=poslin(net.weights{i}*x+net.biases{i});
    end
    y=net.weights{end}*x+net.biases{end};
    
end



