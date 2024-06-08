clear all
clc
close all


load('networks.mat')
model.type='FFNN';
for i=1:400
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

str1=  or_operation( or_operation(  str13 ,str12  ) ,str11  );
str01= G_operation( str1, 1,50);

str21= str2stl('[4,0]');
str22= str2stl('[5,0]');
str23= str2stl('[6,0]');

str2=  or_operation( or_operation(  str23 ,str22 ) ,str21  );
str02= G_operation( str2, 1,50);



str31= str2stl('[7,0]');
str32= str2stl('[8,0]');
str33= str2stl('[9,0]');
str34= str2stl('[10,0]');
str35= str2stl('[11,0]');
str36= str2stl('[12,0]');

str3= and_operation(  and_operation( and_operation( and_operation( and_operation(  str36 ,str35 ) ,str34  ) ,str33 ) ,str32 ) ,str31 );
str03= F_operation( str3, 35,50);


str= and_operation( and_operation(  str03 , str02) , str01);

for i=1:100
    model(i).C=eye(3);
    model(i).D=zeros(3,1);
end


predicate_transformation.C(1,:)=[ -7.15   -14.3  0];
predicate_transformation.d(1,1)=  -1;
predicate_transformation.C(2,:)=[  1   0  0];
predicate_transformation.d(2,1)= -0.3;
predicate_transformation.C(3,:)=[  0  1  0];
predicate_transformation.d(3,1)= 0.12;


predicate_transformation.C(4,:)=[  7.15  9.5  0];
predicate_transformation.d(4,1)=   1;
predicate_transformation.C(5,:)=[  -1   0  0];
predicate_transformation.d(5,1)= 0.1;
predicate_transformation.C(6,:)=[  0  -1  0];
predicate_transformation.d(6,1)= -0.33;



predicate_transformation.C(7,:)=[  1  0  0];
predicate_transformation.d(7,1)= -0.02;
predicate_transformation.C(8,:)=[ -1  0  0];
predicate_transformation.d(8,1)=  0.04;
predicate_transformation.C(9,:)=[  0  1  0];
predicate_transformation.d(9,1)=  0.02;
predicate_transformation.C(10,:)=[ 0 -1  0];
predicate_transformation.d(10,1)= 0.02;
predicate_transformation.C(11,:)=[ 0  0  1];
predicate_transformation.d(11,1)= 0.02;
predicate_transformation.C(12,:)=[ 0  0 -1];
predicate_transformation.d(12,1)= 0.02;

% [Net, border]  = Trapezius_maker(model, nncontroller, predicate_transformation , str);


Center= [0.375 ; -0.325 ; 0.375 ];  %%% This is the center of the set of initial states $\mathcal{X}_0$. 
epsilon=[0.025;   0.025 ; 0.025 ];    %%% This is infinite norm radious of the set of initial states.

tic
num_cores=1;
interval = Trapezius_ReLUplex_verifier(Center, epsilon, model, nncontroller, predicate_transformation, str, analysis_type, num_cores);
Computation_time=toc