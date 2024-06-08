clear all
clc
close all

load('networks.mat')

for i=1:55
    model(i).weights = Model_nn.weights;
    model(i).biases = Model_nn.biases;
    model(i).type = 'FFNN';
    model(i).C = eye(6);
    model(i).D = zeros(6,1);
    model(i).layers = Model_nn.layers;
    
    nncontroller(i) = nn(i);
end


cell_P1_11= str2stl('[1,0]');
cell_P1_12= str2stl('[2,0]');
cell_P1_13= str2stl('[3,0]');
cell_P1_14= str2stl('[4,0]');
cell_P1_15= str2stl('[5,0]');
cell_P1_16= str2stl('[6,0]');


cell_P1 = and_operation(and_operation(and_operation( and_operation(  and_operation( cell_P1_16 , cell_P1_15) , cell_P1_14 ) ,  cell_P1_13 ) , cell_P1_12 ) ,  cell_P1_11 );


cell_P2_11= str2stl('[7,0]');
cell_P2_12= str2stl('[8,0]');
cell_P2_13= str2stl('[9,0]');
cell_P2_14= str2stl('[10,0]');
cell_P2_15= str2stl('[11,0]');
cell_P2_16= str2stl('[12,0]');

cell_P2 = and_operation(and_operation(and_operation( and_operation(  and_operation( cell_P2_16 , cell_P2_15) , cell_P2_14 ) ,  cell_P2_13 ) , cell_P2_12 ) ,  cell_P2_11 );


cell_P3_11= str2stl('[13,0]');
cell_P3_12= str2stl('[14,0]');
cell_P3_13= str2stl('[15,0]');
cell_P3_14= str2stl('[16,0]');
cell_P3_15= str2stl('[17,0]');
cell_P3_16= str2stl('[18,0]');

cell_P3 = and_operation(and_operation(and_operation( and_operation(  and_operation( cell_P3_16 , cell_P3_15) , cell_P3_14 ) ,  cell_P3_13 ) , cell_P3_12 ) ,  cell_P3_11 );




cell_P4_11= str2stl('[19,0]');
cell_P4_12= str2stl('[20,0]');
cell_P4_13= str2stl('[21,0]');
cell_P4_14= str2stl('[22,0]');
cell_P4_15= str2stl('[23,0]');
cell_P4_16= str2stl('[24,0]');

cell_P4 =  or_operation( or_operation( or_operation(  or_operation(  or_operation(  cell_P4_16 , cell_P4_15) , cell_P4_14 ) ,  cell_P4_13 ) , cell_P4_12 ) ,  cell_P4_11 );


% cell_F3 = F_operation(cell_P3,25,50);
cell_1 = or_operation(  cell_P1, cell_P2 );
% cell_0 = and_operation( cell_1 , cell_F3 );

cell_F  = F_operation(cell_1, 1,20);
cell_G  = G_operation(cell_P4,1,20);

cell_STL= and_operation(cell_F , cell_G);


%%% P1
predicate_transformation.C(1,:)=[ -1  0   0  0  0  0 ];
predicate_transformation.d(1,1)= 0.1381;

predicate_transformation.C(2,:)=[ 0  -1   0  0  0  0  ];
predicate_transformation.d(2,1)= 0.3931;

predicate_transformation.C(3,:)=[ 0  0   -1  0  0  0  ];
predicate_transformation.d(3,1)= 0.0531;

predicate_transformation.C(4,:)=[1  0   0  0  0  0  ];
predicate_transformation.d(4,1)= -0.0319;

predicate_transformation.C(5,:)=[ 0 1   0  0  0  0  ];
predicate_transformation.d(5,1)= -0.2869;

predicate_transformation.C(6,:)=[ 0  0  1  0  0  0  ];
predicate_transformation.d(6,1)= 0.0531;



%%% P2
predicate_transformation.C(7,:)=[-1  0  0  0  0  0  ];
predicate_transformation.d(7,1)= 0.3931;

predicate_transformation.C(8,:)=[ 0 -1  0  0  0  0  ];
predicate_transformation.d(8,1)= 0.1381;

predicate_transformation.C(9,:)=[ 0  0 -1  0  0  0  ];
predicate_transformation.d(9,1)= 0.0531;

predicate_transformation.C(10,:)=[ 1  0  0  0  0  0  ];
predicate_transformation.d(10,1)= -0.2869;

predicate_transformation.C(11,:)=[ 0  1  0  0  0  0  ];
predicate_transformation.d(11,1)= -0.0319;

predicate_transformation.C(12,:)=[ 0   0  1  0  0  0  ];
predicate_transformation.d(12,1)= 0.0513;




%%% P3
predicate_transformation.C(13,:)=[-1  0  0  0  0  0  ];
predicate_transformation.d(13,1)= 0.3931;

predicate_transformation.C(14,:)=[ 0 -1  0  0  0  0  ];
predicate_transformation.d(14,1)= 0.3931;

predicate_transformation.C(15,:)=[ 0  0 -1  0  0  0  ];
predicate_transformation.d(15,1)= -0.0744;

predicate_transformation.C(16,:)=[ 1  0  0  0  0  0  ];
predicate_transformation.d(16,1)= -0.2869;

predicate_transformation.C(17,:)=[ 0  1  0  0  0  0  ];
predicate_transformation.d(17,1)= -0.2869;

predicate_transformation.C(18,:)=[ 0   0  1  0  0  0  ];
predicate_transformation.d(18,1)= 0.1806;



%%% P4
predicate_transformation.C(19,:)=[ 1  0  0  0  0  0  ];
predicate_transformation.d(19,1)= -0.2656;

predicate_transformation.C(20,:)=[ 0  1  0  0  0  0  ];
predicate_transformation.d(20,1)= -0.2656;

predicate_transformation.C(21,:)=[ 0  0  1  0  0  0  ];
predicate_transformation.d(21,1)= -0.0531;

predicate_transformation.C(22,:)=[ -1  0  0  0  0  0  ];
predicate_transformation.d(22,1)= 0.1594;

predicate_transformation.C(23,:)=[ 0  -1  0  0  0  0  ];
predicate_transformation.d(23,1)= 0.1594;

predicate_transformation.C(24,:)=[ 0   0 -1  0  0  0  ];
predicate_transformation.d(24,1)= -0.0531;




[net, border]  = Trapezius_maker(model, nncontroller, predicate_transformation , cell_STL);

% load('trapezium.mat')
epsii=0;
lb=[0.0638; 0.0638;  -0.0231; -epsii;  -epsii;  -epsii];
ub=[0.1063; 0.1063 ;  0.0231;  epsii;   epsii;   epsii];
pos{1}=linspace(lb(1),ub(1),8);
pos{2}=linspace(lb(2),ub(2),8);
pos{3}=linspace(lb(3),ub(3),8);

parfor i1=1:8^3
    jj = i1;
    ii1 = mod(jj,8);
    jj = (jj-ii1)/8;
    ii2 = mod( jj , 8);
    jj = (jj-ii2)/8;
    ii3 = mod(jj,8);
    
end
tic
M=2*ones(1,6);
method_Trajectory='Linear-CROWN';
% method_Trajectory='Quadratic-CROWN';
method_STL       ='approx-star' ;
% method_STL       ='exact-star' ;
% LipSDP_type= 'LipSDP-network';
% LipSDP_type= 'LipSDP-layer';
LipSDP_type= 'LipSDP-neuron';
model_dim=[9,10,10,6];
cntrl_dim=[6,10,3,3];
loc_lip  =  Lip_Verify_threshold( net, model_dim, cntrl_dim, border, lb, ub, method_Trajectory, method_STL, LipSDP_type , 0.5 );
[decision, Letsgo, counter_example]  =  Lip_Verify_lip_fixed( net, border, lb, ub, M, method_Trajectory, method_STL, LipSDP_type , loc_lip );
Running_time=toc;