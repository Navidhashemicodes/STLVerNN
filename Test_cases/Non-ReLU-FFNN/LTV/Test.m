clear all
clc
close all
model.type='LTV';
%%%%% be careful about the number of cells per parameter. It should be
%%%%% consistant with the STL specification you will introduce later,


load('LTV_disc.mat')

num_steps=50;

for i=1:num_steps+10
    j=mod(i,30);
    if j==0
        model.A{i}=A_d{30};
        model.B{i}=B_d{30};
    else
        model.A{i}=A_d{j};
        model.B{i}=B_d{j};
    end
    model.E{i}=zeros(2,1);
    model.C{i}=eye(2);
    model.D{i}=zeros(2,1);
    model.type='LTV';
    
end

load('neuro_controller.mat')
for i=1:num_steps+10
    nncontroller(i).weights=net.weights;
    nncontroller(i).biases =net.biases;
    func=cell(7,1);
    func(:)={'tanh'};
    nncontroller(i).layers{1}=func;
    nncontroller(i).layers{2}=func;
end

cell_P1 = str2stl('[1,0]');


cell_P2 = str2stl('[2,0]');


cell_P3 = str2stl('[3,0]');


cell_P4 = str2stl('[4,0]');


cell_P5 = str2stl('[5,0]');


cell_F1 = F_operation(cell_P1,30,35);

cell_F2 = F_operation(cell_P2,37,43);

cell_F3 = F_operation(cell_P3,45,50);

cell_F4 = F_operation(cell_P4,32,38);

cell_G  = G_operation(cell_P5,1 ,50);

cell_STL= and_operation(and_operation(and_operation(and_operation( cell_F1 , cell_F2  ) , cell_F3 ) , cell_F4 ) , cell_G );


predicate_transformation.C(1,:)=[-1 0];
predicate_transformation.d(1,1)=-0.5;
predicate_transformation.C(2,:)=[1 0 ];
predicate_transformation.d(2,1)=0.4;
predicate_transformation.C(3,:)=[-1 0];
predicate_transformation.d(3,1)=0;
predicate_transformation.C(4,:)=[0 1 ];
predicate_transformation.d(4,1)=-1;
predicate_transformation.C(5,:)=[1 -1 ];
predicate_transformation.d(5,1)=5.5;



[net, border]  = Trapezius_maker(model, nncontroller, predicate_transformation , cell_STL);

lb=[-1;-1];
ub=[ 1; 1];
% 
% NNN=30;
% 
% xs = linspace(lb(1),ub(1),NNN);
% ys = linspace(lb(2),ub(2),NNN);
% 
% 
% robs = zeros(NNN,NNN);
% for xx = 1:length(xs)
%     for yy = 1:length(ys)
%         robs(xx,yy) = NN(net,[xs(xx);ys(yy)]);
%     end
% end
% 
% figure;
% view(3);
% surface(xs,ys,robs);





method_Trajectory='Linear-CROWN';
% method_Trajectory='Quadratic-CROWN';
method_STL       ='approx-star' ;
% method_STL       ='exact-star' ;
LipSDP_type= 'LipSDP-layer';

loc_lip=zeros(H,H);
status=zeros(H,H);
Run_time1=zeros(H,H);
decision=zeros(H,H);
Letsgo=zeros(H,H);
counter_example=zeros(H,H,2);
Run_time2=zeros(H,H);

parfor i=1:H
    for j=1:H
        tic
        [loc_lip(i,j), status(i,j)]  =  Lip_Verify_threshold( net, [3,2], [2,7,7,1], border, lb, ub, method_Trajectory, method_STL, LipSDP_type , [1,1] );
        Run_time1(i,j)=toc;
        tic
        [decision(i,j), Letsgo(i,j), counter_example(i,j,:)]  =  Lip_Verify_lip_fixed( net,  border, lb, ub, [2 , 2], method_Trajectory, method_STL, LipSDP_type , loc_lip );
        Run_time2(i,j)=toc;
        save(['result_for' num2str(i) 'and' num2str(j)] , 'Run_time1', 'Run_time2', 'loc_lip', 'status', 'decision', 'counter_example')
    end
end