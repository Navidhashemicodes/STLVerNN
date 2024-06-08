clear all
clc
close all

load('Toy_example_2.mat')

nncontroller.weights=NN_actor.weights;
model.weights=NN_model.weights;

nncontroller.biases=NN_actor.Biases;
model.biases=NN_model.Biases;


L=cell(5,1);
L(:)={'tanh'};
nncontroller.layers{1}=L;

L=cell(5,1);
L(:)={'tanh'};
model.layers{1}=L;

model.type='FFNN';

for i=1:50
    model.C{i}=eye(2);
    model.D{i}=zeros(2,1);
end


cell_P1_11= str2stl('[1,0]');
cell_P1_12= str2stl('[2,0]');
cell_P1_13= str2stl('[3,0]');
cell_P1_14= str2stl('[4,0]');

cell_P1 = and_operation( and_operation(  and_operation( cell_P1_14 , cell_P1_13) , cell_P1_12 ) ,  cell_P1_11 );


cell_P2_11= str2stl('[5,0]');
cell_P2_12= str2stl('[6,0]');
cell_P2_13= str2stl('[7,0]');
cell_P2_14= str2stl('[8,0]');

cell_P2 = and_operation( and_operation(  and_operation( cell_P2_14 , cell_P2_13) , cell_P2_12 ) ,  cell_P2_11 );

cell_F1 = F_operation(cell_P2,9,13);

cell_1 = and_operation(cell_P1, cell_F1);
cell_STL= F_operation(cell_1 , 3,6);



%%% P1
x0=-0.4;y0=-0.7;x1= 0.3;y1=-0.4;
a=-(y1-y0)/(x1-x0); b=(x0*(y1-y0)/(x1-x0))-y0;
predicate_transformation.C(1,:)=[ a  1 ];
predicate_transformation.d(1,1)= b;


x0=-0.8;y0=-0.5;x1= -0.4;y1=-0.7;
a=-(y1-y0)/(x1-x0); b=(x0*(y1-y0)/(x1-x0))-y0;
predicate_transformation.C(2,:)=[ a  1 ];
predicate_transformation.d(2,1)= b;


x0=-0.8;y0=-0.5;x1= -0.1;y1=-0.2;
a=-(y1-y0)/(x1-x0); b=(x0*(y1-y0)/(x1-x0))-y0;
predicate_transformation.C(3,:)=[ -a  -1 ];
predicate_transformation.d(3,1)= -b;


x0=-0.1;y0=-0.2;x1= 0.3;y1=-0.4;
a=-(y1-y0)/(x1-x0); b=(x0*(y1-y0)/(x1-x0))-y0;
predicate_transformation.C(4,:)=[ -a  -1 ];
predicate_transformation.d(4,1)= -b;


%%% P2
predicate_transformation.C(5,:)=[ 1  0 ];
predicate_transformation.d(5,1)= 1.1;

predicate_transformation.C(6,:)=[-1  0 ];
predicate_transformation.d(6,1)= -0.7;

predicate_transformation.C(7,:)=[ 0 -1 ];
predicate_transformation.d(7,1)= 0;

predicate_transformation.C(8,:)=[ 0  1 ];
predicate_transformation.d(8,1)= 0.3;

[net, border]  = Trapezius_maker(model, nncontroller, predicate_transformation , cell_STL);
M=2;
lb=[1;1];
ub=[2;2];
tic
% [decision, Letsgo]= Lip_Verify2(net, border, lb,ub, M );

method_Trajectory='Linear-CROWN';
method_STL       ='approx-star' ;
LipSDP_type= 'LipSDP-layer';
[decision, Letsgo, counter_example]  =  Lip_Verify( net, border, lb, ub, [2 , 2], method_Trajectory, method_STL, LipSDP_type );
Running_time=toc;
