clear all
clc
close all


load('final_example_BigB.mat')
load('network.mat')



model.type='LTV';

cell_P1_1=str2stl('[1,0]');
cell_P1_2=str2stl('[2,0]');
cell_P1_3=str2stl('[3,0]');
cell_P1_4=str2stl('[4,0]');

cell_P1_g = and_operation( cell_P1_4 , and_operation( cell_P1_3 , and_operation( cell_P1_2 , cell_P1_1 ) ) );
cell_P1   = G_operation(   cell_P1_g , 9,16);

cell_P2_1=str2stl('[5,0]');
cell_P2_2=str2stl('[6,0]');
cell_P2_3=str2stl('[7,0]');
cell_P2_4=str2stl('[8,0]');

cell_P2_g = and_operation( cell_P2_4 , and_operation( cell_P2_3 , and_operation( cell_P2_2 , cell_P2_1 ) ) );
cell_P2   = G_operation(   cell_P2_g , 17,24);
cell_STL  = and_operation( cell_P1 , cell_P2 );


cell_P3_1=str2stl('[9,0]');
cell_P3_2=str2stl('[10,0]');
cell_P3_3=str2stl('[11,0]');
cell_P3_4=str2stl('[12,0]');

cell_P3_g = and_operation( cell_P3_4 , and_operation( cell_P3_3 , and_operation( cell_P3_2 , cell_P3_1 ) ) );
cell_P3   = G_operation(   cell_P3_g , 25,32);
cell_STL  = and_operation( cell_STL , cell_P3 );



cell_P4_1=str2stl('[13,0]');
cell_P4_2=str2stl('[14,0]');
cell_P4_3=str2stl('[15,0]');
cell_P4_4=str2stl('[16,0]');

cell_P4_g = and_operation( cell_P4_4 , and_operation( cell_P4_3 , and_operation( cell_P4_2 , cell_P4_1 ) ) );
cell_P4   = G_operation(   cell_P4_g , 33,40);
cell_STL  = and_operation( cell_STL , cell_P4 );


cell_P5_1=str2stl('[17,0]');
cell_P5_2=str2stl('[18,0]');
cell_P5_3=str2stl('[19,0]');
cell_P5_4=str2stl('[20,0]');

cell_P5_g = and_operation( cell_P5_4 , and_operation( cell_P5_3 , and_operation( cell_P5_2 , cell_P5_1 ) ) );
cell_P5   = G_operation(   cell_P5_g , 41,43);
cell_STL  = and_operation( cell_STL , cell_P5 );

cell_P6_1=str2stl('[21,0]');
cell_P6_2=str2stl('[22,0]');
cell_P6_3=str2stl('[23,0]');
cell_P6_4=str2stl('[24,0]');

cell_P6_g = and_operation( cell_P6_4 , and_operation( cell_P6_3 , and_operation( cell_P6_2 , cell_P6_1 ) ) );
cell_P6   = G_operation(   cell_P6_g , 44,60);
cell_STL  = and_operation( cell_STL , cell_P6 );



model.type='LTV';

for i=1:200
    model.A{i}=A;
    model.B{i}=B;
    model.E{i}=zeros(2,1);
    model.C{i}=eye(2);
    model.D{i}=zeros(2,1);
    

    NN_controller(i).weights=nn_controller.weights;
    NN_controller(i).biases =nn_controller.biases;
    NN_controller(i).layers =nn_controller.layers;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%  P1
b1=-250;a1=35;b2=180;a2=-45;
% y=x/3+a1;
predicate_transformation.C(1,:)=[1/3 -1];
predicate_transformation.d(1,1)=a1;

% y=x/3+a2;
predicate_transformation.C(2,:)=[-1/3 1];
predicate_transformation.d(2,1)=-a2;

% y=-3*x+b1;
predicate_transformation.C(3,:)=[3 1];
predicate_transformation.d(3,1)=-b1;

% y=-3*x+b2;
predicate_transformation.C(4,:)=[-3 -1];
predicate_transformation.d(4,1)=b2;



%%%%  P2
b1=-130;a1=16;b2=80;a2=-22;
% y=x/3+a1;
predicate_transformation.C(5,:)=[1/3 -1];
predicate_transformation.d(5,1)=a1;

% y=x/3+a2;
predicate_transformation.C(6,:)=[-1/3 1];
predicate_transformation.d(6,1)=-a2;

% y=-3*x+b1;
predicate_transformation.C(7,:)=[3 1];
predicate_transformation.d(7,1)=-b1;

% y=-3*x+b2;
predicate_transformation.C(8,:)=[-3 -1];
predicate_transformation.d(8,1)=b2;



%%%%  P3
b1=-55;a1=7;b2=40;a2=-10;
% y=x/3+a1;
predicate_transformation.C(9,:)=[1/3 -1];
predicate_transformation.d(9,1)=a1;

% y=x/3+a2;
predicate_transformation.C(10,:)=[-1/3 1];
predicate_transformation.d(10,1)=-a2;

% y=-3*x+b1;
predicate_transformation.C(11,:)=[3 1];
predicate_transformation.d(11,1)=-b1;

% y=-3*x+b2;
predicate_transformation.C(12,:)=[-3 -1];
predicate_transformation.d(12,1)=b2;



%%%%  P4
b1=-25;a1=3;b2=15;a2=-4;
% y=x/3+a1;
predicate_transformation.C(13,:)=[1/3 -1];
predicate_transformation.d(13,1)=a1;

% y=x/3+a2;
predicate_transformation.C(14,:)=[-1/3 1];
predicate_transformation.d(14,1)=-a2;

% y=-3*x+b1;
predicate_transformation.C(15,:)=[3 1];
predicate_transformation.d(15,1)=-b1;

% y=-3*x+b2;
predicate_transformation.C(16,:)=[-3 -1];
predicate_transformation.d(16,1)=b2;




%%%%  P5
b1=-9;a1=1.1;b2=5;a2=-1.2;
% y=x/3+a1;
predicate_transformation.C(17,:)=[1/3 -1];
predicate_transformation.d(17,1)=a1;

% y=x/3+a2;
predicate_transformation.C(18,:)=[-1/3 1];
predicate_transformation.d(18,1)=-a2;

% y=-3*x+b1;
predicate_transformation.C(19,:)=[3 1];
predicate_transformation.d(19,1)=-b1;

% y=-3*x+b2;
predicate_transformation.C(20,:)=[-3 -1];
predicate_transformation.d(20,1)=b2;



%%%%  P6
b1=-1.2;a1=0.7;b2=3.5;a2=-0.6;
% y=x/3+a1;
predicate_transformation.C(21,:)=[1/3 -1];
predicate_transformation.d(21,1)=a1;

% y=x/3+a2;
predicate_transformation.C(22,:)=[-1/3 1];
predicate_transformation.d(22,1)=-a2;

% y=-3*x+b1;
predicate_transformation.C(23,:)=[3 1];
predicate_transformation.d(23,1)=-b1;

% y=-3*x+b2;
predicate_transformation.C(24,:)=[-3 -1];
predicate_transformation.d(24,1)=b2;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





Center= [-41;88];  %%% This is the center of the set of initial states $\mathcal{X}_0$. 
epsilon=[  1; 1 ];    %%% This is infinite norm radious of the set of initial states.
analysis_type='exact-star';
% analysis_type='approx-star';

num_Core=1;
tic
interval = Trapezius_ReLUplex_verifier(Center, epsilon, model, NN_controller, predicate_transformation, cell_STL, analysis_type, num_Core);
Computation_time=toc;