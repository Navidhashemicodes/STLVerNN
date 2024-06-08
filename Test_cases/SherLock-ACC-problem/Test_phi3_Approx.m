clear all
clc
close all





load('networks.mat')
model.type='FFNN';
for i=1:100
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

% str = G_operation(str2stl('[1,0]'),0,50);
str = G_operation(    or_operation(     str2stl('[1,0]')  ,    F_operation(str2stl('[2,0]'),0,3)    ),         1,50        );

E=[ 0  0  0  0  1  0;...
    1  0  0 -1  0  0;...
    0  1  0  0 -1  0];
t_gap=1.4;
V_set=30;
for i=1:100
    model(i).C=[zeros(2,6) ; E ];
    model(i).D=[V_set ; t_gap; 0; 0; 0 ];
end

D_default=10;
D_default2=12;

predicate_transformation.C(1,:)=[1 0 0 -1 -t_gap 0 ];
predicate_transformation.d(1,1)=-D_default;
predicate_transformation.C(2,:)=[1 0 0 -1 -t_gap 0 ];
predicate_transformation.d(2,1)=-D_default2;

Center=[100;32.1;0;10.5;30.1;0];  %%% This is the center of the set of initial states $\mathcal{X}_0$. 
epsilon=[10;0.1;0;0.5;0.1;0];    %%% This is infinite norm radious of the set of initial states.

num_Core=1;
tic
interval = Trapezius_ReLUplex_verifier(Center, epsilon, model, nncontroller, predicate_transformation, str, analysis_type, num_Core);
Computation_time=toc