clear all
clc
close all


A_c=[ 0  1  0  0  0  0  ;...
      0  0  1  0  0  0  ;...
      0  0 -2  0  0  0  ;...
      0  0  0  0  1  0  ;...
      0  0  0  0  0  1  ;...
      0  0  0  0  0 -2  ];
B_c=[ 0  0  ;...
      0  0  ;...
      2  0  ;...
      0  0  ;...
      0  0  ;...
      0  2  ];
  
V_set=30;
t_step=0.1;
t_gap=1.4;
D_default=10;
D_default2=12;

E_tran=[ 0  0  0  0  1  0;...
         1  0  0 -1  0  0;...
         0  1  0  0 -1  0];

A= expm(A_c*t_step);
B= G_maker(t_step,A_c,B_c);



model.type='LTV';
%%%%% be careful about the number of cells per parameter. It should be
%%%%% consistant with the STL specification you will introduce later,


% str=input('please insert the STL specification based on manual: ');
% str = G_operation(F_operation(str2stl('[1,0]'),0,3),1,48);
str = G_operation(    or_operation(     str2stl('[1,0]')  ,    F_operation(str2stl('[2,0]'),0,3)    ),         1,50        );

load('nn3_20.mat')
for i=1:100
    model.A{i}=A;
    model.B{i}=B(:,2);
    model.E{i}=B(:,1)*(-2);
    model.C{i}=[zeros(2,6) ; E_tran ];
    model.D{i}=[V_set ; t_gap; 0; 0; 0 ];
    model.type='LTV';
    
    nncontroller(i).weights=network.weights;
    nncontroller(i).biases =network.bias;
    func=cell(20,1);
    func(:)={'poslin'};
    nncontroller(i).layers{1}=func;
    nncontroller(i).layers{2}=func;
    nncontroller(i).layers{3}=func;
end




predicate_transformation.C(1,:)=[1 0 0 -1 -t_gap 0 ];
predicate_transformation.d(1,1)=-D_default;
predicate_transformation.C(2,:)=[1 0 0 -1 -t_gap 0 ];
predicate_transformation.d(2,1)=-D_default2;

Center=[100;32.1;0;10.5;30.1;0];  %%% This is the center of the set of initial states $\mathcal{X}_0$. 
epsilon=[10;0.1;0;0.5;0.1;0];    %%% This is infinite norm radious of the set of initial states.
% analysis_type='exact-star';
analysis_type='approx-star';

num_Core=1;
tic
interval = Trapezius_ReLUplex_verifier(Center, epsilon, model, nncontroller, predicate_transformation, str, analysis_type, num_Core);
Computation_time=toc