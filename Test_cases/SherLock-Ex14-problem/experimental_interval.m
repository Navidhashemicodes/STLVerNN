clear all
clc
close all

load('networks.mat')
model=Model_nn;
nncontroller=controller_nn;
model.type='FFNN';
%%%%% be careful about the number of cells per parameter. It should be
%%%%% consistant with the STL specification you will introduce later,


analysis_type='exact-star';
% analysis_type='approx-star';


str11= str2stl('[1,0]');
str12= str2stl('[2,0]');
str13= str2stl('[3,0]');
str14= str2stl('[4,0]');

str1 = and_operation( and_operation(  and_operation( str13 , str14) , str12 ) ,  str11 );
str01= F_operation(str1, 75, 100);

p_1 = 'p1 /\ p2 /\ p3 /\ p4';

str21= str2stl('[5,0]');
str22= str2stl('[6,0]');
str23= str2stl('[7,0]');

str2= or_operation(  or_operation(  str23, str22 ) , str21  );
str02= G_operation( str2, 1,100);

p_2 = 'p5 \/ p6 \/ p7 ';

str31= str2stl('[8,0]');
str32= str2stl('[9,0]');
str33= str2stl('[10,0]');

str3= or_operation(  or_operation(  str33, str32 ) , str31  );
str03= G_operation( str3, 1,100);

p_3 ='p8 \/ p9  \/  p10 ';


str= and_operation( and_operation( str03, str02) , str01);

p = '( <>_[75,101] ( p1 /\ p2 /\ p3 /\ p4 ) )  /\  ( []_[2,101] ( p5 \/ p6 \/ p7 ) )  /\  ( []_[2,101] ( p8 \/ p9 \/ p10 ) )';




for i=1:200
    model.C{i}=eye(2);
    model.D{i}=zeros(2,1);
end


predicate_transformation.C(1,:)=[ 1  0 ];
predicate_transformation.d(1,1)= -0.1;
predicate_transformation.C(2,:)=[-1  0 ];
predicate_transformation.d(2,1)= 0.15;
predicate_transformation.C(3,:)=[ 0 -1 ];
predicate_transformation.d(3,1)= 0.02;
predicate_transformation.C(4,:)=[ 0  1 ];
predicate_transformation.d(4,1)= 0.02;


p_map(1).str= 'p1';
p_map(1).A  = [-1 0];
p_map(1).b  =  -0.1;

p_map(2).str= 'p2';
p_map(2).A  = [ 1 0];
p_map(2).b  =  0.15;

p_map(3).str= 'p3';
p_map(3).A  = [ 0 1];
p_map(3).b  =  0.2;

p_map(4).str= 'p4';
p_map(4).A  = [ 0 -1];
p_map(4).b  =  0.2;

predicate_transformation.C(5,:)=[ 2  -1 ];
predicate_transformation.d(5,1)= -1;
predicate_transformation.C(6,:)=[-1  0 ];
predicate_transformation.d(6,1)= 0.5;
predicate_transformation.C(7,:)=[ 0  1 ];
predicate_transformation.d(7,1)= -0.3;

p_map(5).str= 'p5';
p_map(5).A  = [ -2 1];
p_map(5).b  =  -1;

p_map(6).str= 'p6';
p_map(6).A  = [ 1 0];
p_map(6).b  =  0.5;


p_map(7).str= 'p7';
p_map(7).A  = [ 0 -1];
p_map(7).b  =  -0.3;



predicate_transformation.C(8,:)=[ -3  3 ];
predicate_transformation.d(8,1)= 2.1;
predicate_transformation.C(9,:)=[ 1  0 ];
predicate_transformation.d(9,1)= -0.8;
predicate_transformation.C(10,:)=[ 0 -1 ];
predicate_transformation.d(10,1)= 0;


p_map(8).str= 'p8';
p_map(8).A  = [ 3 -3];
p_map(8).b  =  2.1;

p_map(9).str= 'p9';
p_map(9).A  = [ -1 0];
p_map(9).b  =  -0.8;


p_map(10).str= 'p10';
p_map(10).A  = [ 0 1];
p_map(10).b  =  0;





Center=[0.85;0.45];  %%% This is the center of the set of initial states $\mathcal{X}_0$. 
epsilon=[0.05;0.05];    %%% This is infinite norm radious of the set of initial states.

[net, border] = Trapezius_maker(model, nncontroller, predicate_transformation , str);
Net=Linear_Nonlinear(net, border);


for kk=1:40
s(:,1)=Center+(2*rand(2,1)-1).*epsilon;
for i=1:100
    a(:,i)=pred(nncontroller,s(:,i));
    s(:,i+1)=pred(model, [s(:,i);a(:,i)]);
end

seqT=1:101;
seqS=s';
[rob,~] = dp_taliro(p,p_map,seqS,seqT');


G1=max((predicate_transformation.C(5:7,:)*s(:,1)+predicate_transformation.d(5:7,1))');
G2=max((predicate_transformation.C(8:10,:)*s(:,1)+predicate_transformation.d(8:10,1))');

for i=2:101
    G1=min(G1, max((predicate_transformation.C(5:7,:)*s(:,i)+predicate_transformation.d(5:7,1)))');
    G2=min(G2, max((predicate_transformation.C(8:10,:)*s(:,i)+predicate_transformation.d(8:10,1)))');
end
F1=min((predicate_transformation.C(1:4,:)*s(:,1)+predicate_transformation.d(1:4,1))');
for i=75:101
    F1=max(F1, min((predicate_transformation.C(1:4,:)*s(:,i)+predicate_transformation.d(1:4,1)))');
end
clc
rho1(kk)=min([F1,G1,G2]);
rho2=NN(Net, s(:,1));
disp(rob)
disp(rho1(end))
disp(rho2)
pause(2)
clc
end 
% 
% rho2=NN(Net, s(:,1))
% rho3=NN(net, s(:,1))
% figure
% plot(rho1)   
    
function y = pred(net, x)
   
    len=length(net.weights)-1;
    for i=1:len
        x=poslin(net.weights{i}*x+net.biases{i});
    end
    y=net.weights{end}*x+net.biases{end};
    
end