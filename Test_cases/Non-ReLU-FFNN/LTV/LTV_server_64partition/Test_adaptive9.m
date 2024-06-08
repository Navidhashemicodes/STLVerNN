clear all
clc
close all
model.type='LTV';

%%%%% be careful about the number of cells per parameter. It should be
%%%%% consistant with the STL specification you will introduce later,

addpath(genpath('/home1/navidhas/Personal_Trapezius_toolbox_time_variant'));
addpath(genpath('/home1/navidhas/nnv'));
addpath(genpath('/home1/navidhas/Mosek_Home'));
addpath(genpath('/home1/navidhas/mosek'));


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

% lb=[-1   ; -1  ];
% ub=[0 ; 0];

Lb=[-1; -1];
Ub=[ 0; 0 ];
n1=8;
Points=[linspace(Lb(1), Ub(1) , n1+1); linspace(Lb(2), Ub(2) , n1+1)];
epsi=(Ub-Lb)/n1;
init_dim=2;

hh=9;

hh=hh-1;
rr=dec2base(hh,n1);
len=length(rr);
R=zeros(1,init_dim);
if len==1
    R(1)=0;
    R(2)=0;
    R(3)=str2num(rr(1));
elseif len==2
    R(1)=0;
    R(2)=str2num(rr(1));
    R(3)=str2num(rr(2));
elseif len==3
    R(1)=str2num(rr(1));
    R(2)=str2num(rr(2));
    R(3)=str2num(rr(3));
end
v=R+1;
lb=[Points(1,v(1)); Points(2,v(2))];
ub=lb+epsi;




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
M=[2,2];
tic
[Lb, Ub, alpha_params, beta_params] =  Trapezius_bounds_General( 0.5*(lb+ub), 0.5*(ub-lb), net, border,  method_Trajectory, method_STL);
[Lip, status, time]  =  Trapezius_Lip_SDP( net, border, alpha_params, beta_params , Lb, Ub,  LipSDP_type );
% status=1;
% Lip=[];
[decision, Letsgo, counter_example]  =  Lip_Verify_adaptive1( net, border, lb, ub, M, method_Trajectory, method_STL, LipSDP_type, Lip , status  );
run_time=toc;

save(['Results_hh_' num2str(hh+1)], "run_time", "decision","counter_example", "lb", "ub")