clear all
clc
close all

load('Best_example.mat')
model.type='FFNN';
for i=1:50
    nncontroller(i).weights=NN_actor.weights;
    model(i).weights=NN_model.weights;
    
    nncontroller(i).biases=NN_actor.Biases;
    model(i).biases=NN_model.Biases;
    
    
    L=cell(8,1);
    L(:)={'poslin'};
    nncontroller(i).layers{1}=L;
    
    L=cell(10,1);
    L(:)={'poslin'};
    model(i).layers{1}=L;
  
end
%%%%% be careful about the number of cells per parameter. It should be
%%%%% consistant with the STL specification you will introduce later,


analysis_type='exact-star';
% analysis_type='approx-star';

cell_P1_11= str2stl('[1,0]');
cell_P1_12= str2stl('[2,0]');
cell_P1_13= str2stl('[3,0]');
cell_P1_14= str2stl('[4,0]');

cell_P1 = and_operation( and_operation(  and_operation( cell_P1_13 , cell_P1_14) , cell_P1_12 ) ,  cell_P1_11 );


cell_P2_11= str2stl('[5,0]');
cell_P2_12= str2stl('[6,0]');
cell_P2_13= str2stl('[7,0]');
cell_P2_14= str2stl('[8,0]');

cell_P2 = and_operation( and_operation(  and_operation( cell_P2_13 , cell_P2_14) , cell_P2_12 ) ,  cell_P2_11 );




cell_P3_11= str2stl('[9,0]');
cell_P3_12= str2stl('[10,0]');
cell_P3_13= str2stl('[11,0]');
cell_P3_14= str2stl('[12,0]');

cell_P3 = and_operation( and_operation(  and_operation( cell_P3_13 , cell_P3_14) , cell_P3_12 ) ,  cell_P3_11 );





cell_P4_11= str2stl('[13,0]');
cell_P4_12= str2stl('[14,0]');
cell_P4_13= str2stl('[15,0]');
cell_P4_14= str2stl('[16,0]');

cell_P4 = and_operation( and_operation(  and_operation( cell_P4_13 , cell_P4_14) , cell_P4_12 ) ,  cell_P4_11 );



% %%% Reject it
% cell_1st = F_operation(cell_P4,    9,10);
% cell_2nd = and_operation(cell_1st, cell_P3);
% cell_3rd = F_operation(cell_2nd,    6,7);
% cell_4th = and_operation(cell_3rd, cell_P2);
% cell_5th = F_operation(cell_4th,    6, 11);
% cell_6th = and_operation(cell_5th, cell_P1);
% cell_7th = F_operation(cell_6th,    5, 8);
% cell_STL = cell_7th;


%% vaerify it
cell_1st = F_operation(cell_P4,    8,9);
cell_2nd = and_operation(cell_1st, cell_P3);
cell_3rd = F_operation(cell_2nd,    6,7);
cell_4th = and_operation(cell_3rd, cell_P2);
cell_5th = F_operation(cell_4th,    6, 11);
cell_6th = and_operation(cell_5th, cell_P1);
cell_7th = F_operation(cell_6th,    5, 8);
cell_STL = cell_7th;

for i=1:200
    model(i).C=eye(2);
    model(i).D=zeros(2,1);
end


%%%P1
predicate_transformation.C(1,:)=[ 1  0 ];
predicate_transformation.d(1,1)= 4;
predicate_transformation.C(2,:)=[-1  0 ];
predicate_transformation.d(2,1)= -2;
predicate_transformation.C(3,:)=[ 0 -1 ];
predicate_transformation.d(3,1)= -0.2;
predicate_transformation.C(4,:)=[ 0  1 ];
predicate_transformation.d(4,1)= 0.8;


%%% P2
predicate_transformation.C(5,:)=[ 1  0 ];
predicate_transformation.d(5,1)= 10;
predicate_transformation.C(6,:)=[-1  0 ];
predicate_transformation.d(6,1)= -7;
predicate_transformation.C(7,:)=[ 0 -1 ];
predicate_transformation.d(7,1)= -0.2;
predicate_transformation.C(8,:)=[ 0  1 ];
predicate_transformation.d(8,1)= 0.8;

%%% P3
predicate_transformation.C(9,:)= [ 1  0 ];
predicate_transformation.d(9,1)=  24;
predicate_transformation.C(10,:)=[-1  0 ];
predicate_transformation.d(10,1)= -14;
predicate_transformation.C(11,:)=[ 0  -1 ];
predicate_transformation.d(11,1)= -0.2;
predicate_transformation.C(12,:)=[ 0 1 ];
predicate_transformation.d(12,1)= 0.6;



%%%P4
predicate_transformation.C(13,:)=[ 1  0 ];
predicate_transformation.d(13,1)=  75;
predicate_transformation.C(14,:)=[-1  0 ];
predicate_transformation.d(14,1)=  -60;
predicate_transformation.C(15,:)=[ 0  -1 ];
predicate_transformation.d(15,1)=  0.6;
predicate_transformation.C(16,:)=[ 0   1 ];
predicate_transformation.d(16,1)= -0.2;





Center= [1.5;1.5];  %%% This is the center of the set of initial states $\mathcal{X}_0$. 
epsilon=[0.5;0.5];    %%% This is infinite norm radious of the set of initial states.

[Net, border]  = Trapezius_maker(model, nncontroller, predicate_transformation , cell_STL);
% Net=Linear_Nonlinear(Net, border);

num_Cores=1;
tic
interval = Trapezius_ReLUplex_verifier(Center, epsilon, model, nncontroller, predicate_transformation, cell_STL, analysis_type,num_Cores);
Computation_time=toc



% l0= Center-epsilon;
% u0= Center+epsilon;
% M=1;
% pos1s=linspace(l0(1), u0(1),M+1);
% pos2s=linspace(l0(2), u0(2),M+1);
% 
% the_center= cell(M,M);
% the_eps   = cell(M,M);
% interval=cell(M,M);
% Computation_time=zeros(M,M);
% [Net, border]  = Trapezius_maker(model, nncontroller, predicate_transformation , cell_STL);
% Net=Linear_Nonlinear(Net, border);
% for i=1:5000
%     s0=[1+rand; 1+rand];
%     R(i)=NN(Net,s0);
% end
% min(R)




% ii=1;
% jj=0;
% for i=ii:M  %%% 5-1 is the number of segements
%     for j=1:M
%         if i~=ii || j>jj
%             the_center{i,j}=0.5*([pos1s(i); pos2s(j)] + [pos1s(i+1); pos2s(j+1)]);
%             the_eps{i,j}   =0.5*([pos1s(i+1); pos2s(j+1)] - [pos1s(i); pos2s(j)]);
%             num_Cores=1;
%             tic
%             interval{i,j} = Trapezius_ReLUplex_verifier_withNet(the_center{i,j}, the_eps{i,j}, Net, analysis_type, num_Cores);
%             Computation_time(i,j)=toc;
%             clc
%         end
%     end
% end
