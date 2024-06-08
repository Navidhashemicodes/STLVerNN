clear all
clc
close all

load('networks.mat')
model.type='FFNN';
for i=1:200
    nncontroller(i).weights=controller_nn.weights;
    model(i).weights=Model_nn.weights;
    
    nncontroller(i).biases=controller_nn.biases;
    model(i).biases=Model_nn.biases;
       
    nncontroller(i).layers=controller_nn.layers;
    model(i).layers=Model_nn.layers;
  
end

%%%%% be careful about the number of cells per parameter. It should be
%%%%% consistant with the STL specification you will introduce later,


% analysis_type='exact-star';
analysis_type='approx-star';


str11= str2stl('[1,0]');
str12= str2stl('[2,0]');
str13= str2stl('[3,0]');
str14= str2stl('[4,0]');

str1 = and_operation( and_operation(  and_operation( str13 , str14) , str12 ) ,  str11 );
str01= F_operation(str1, 75, 100);


str21= str2stl('[5,0]');
str22= str2stl('[6,0]');
str23= str2stl('[7,0]');

str2= or_operation(  or_operation(  str23, str22 ) , str21  );
str02= G_operation( str2, 1,100);

str31= str2stl('[8,0]');
str32= str2stl('[9,0]');
str33= str2stl('[10,0]');

str3= or_operation(  or_operation(  str33, str32 ) , str31  );
str03= G_operation( str3, 1,100);

str= and_operation( and_operation( str03, str02) , str01);





for i=1:200
    model(i).C=eye(2);
    model(i).D=zeros(2,1);
end


predicate_transformation.C(1,:)=[ 1  0 ];
predicate_transformation.d(1,1)= -0.1;
predicate_transformation.C(2,:)=[-1  0 ];
predicate_transformation.d(2,1)= 0.15;
predicate_transformation.C(3,:)=[ 0 -1 ];
predicate_transformation.d(3,1)= 0.02;
predicate_transformation.C(4,:)=[ 0  1 ];
predicate_transformation.d(4,1)= 0.02;


predicate_transformation.C(5,:)=[ 2  -1 ];
predicate_transformation.d(5,1)= -1;
predicate_transformation.C(6,:)=[-1  0 ];
predicate_transformation.d(6,1)= 0.5;
predicate_transformation.C(7,:)=[ 0  1 ];
predicate_transformation.d(7,1)= -0.3;


predicate_transformation.C(8,:)=[ -3  3 ];
predicate_transformation.d(8,1)= 2.1;
predicate_transformation.C(9,:)=[ 1  0 ];
predicate_transformation.d(9,1)= -0.8;
predicate_transformation.C(10,:)=[ 0 -1 ];
predicate_transformation.d(10,1)= 0;


Center=[0.85;0.45];  %%% This is the center of the set of initial states $\mathcal{X}_0$. 
epsilon=[0.05;0.05];    %%% This is infinite norm radious of the set of initial states.

[Net, border]  = Trapezius_maker(model, nncontroller, predicate_transformation , str);
% Net=Linear_Nonlinear(Net, border);


tic
num_cores=1;
interval = Trapezius_ReLUplex_verifier(Center, epsilon, model, nncontroller, predicate_transformation, str, analysis_type, num_cores);
Computation_time=toc