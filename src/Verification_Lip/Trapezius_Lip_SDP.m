function [Lip, status, time] = Trapezius_Lip_SDP( Trapezius_nn, border, alpha_params, beta_params , Lb, Ub,  LipSDP_type )

%%%% This function receives the Tapezium-FFNN and its computed pre-activation bounds (Lb,Ub), alpha_params and beta_params are the lower and upper bounds for the slope of activation functions.
%%%% Although alpha_params and beta_params are computable with Lb,Ub but we also introduce them for convenience in computation. The input border is also the first position of the layers for STL-FFNN in Trapezium-FFNN. 

Net = Linear_Nonlinear(Trapezius_nn , border);  %%% plays with weights and biases to separate the location of linear activations and non-linear activations.
system_dim = size(Net.weights{1},2); %%% this is the dimension of the model. or in another word the dimension of the states.

len = length(Net.layers);
LMI_size = system_dim;   %%% The following for loop computes the LMI size. What is the LMI_size? This is the number of all the nonlinear activation functions + the dimension of the model. 
nl_index = zeros(1, len);  %%% this array stores the position of first nonlinear activation on every single layer. 

for i = 1:len
    nl_index(i) = length(Net.layers{i})-length(alpha_params{i}) +1 ;
    X(i) = length(alpha_params{i});
    LMI_size = LMI_size+X(i);
end
%%% We devide each weight matrix into 4 parts. consider layer \ell-1 and layer \ell the weight matrix between the two is the W{\ell}, layer \ell-1  contains a set of linear activations and the remainings are nonlinear 
%%% activations. The same holds for layer \ell. The part of W{\ell} that connects the linear part of layer \ell-1 to layer \ell is W_L_L{i}. and the part which connects nonlinear part of layer \ell-1 to linear part of 
%%% layer \ell is called W_n_L{i}, the parat that connects nonlinears of layer \ell-1 to nonlinears of layer \ell is called W_n_n{i}. The followin for loop computes these partitions. 

W_L_L{1} = [];   
W_n_L{1} = Net.weights{1}(   1:nl_index(1)-1     ,  :    );
W_L_n{1} = [];
W_n_n{1} = Net.weights{1}(   nl_index(1):end     ,  :    );

for i = 2:len
    W_L_L{i} = Net.weights{i}(   1:nl_index(i)-1     ,     1:nl_index(i-1)-1  );
    W_n_L{i} = Net.weights{i}(   1:nl_index(i)-1     ,     nl_index(i-1):end  );
    W_L_n{i} = Net.weights{i}(   nl_index(i):end     ,     1:nl_index(i-1)-1  );
    W_n_n{i} = Net.weights{i}(   nl_index(i):end     ,     nl_index(i-1):end  );
end


%%% The matrix E_Q{\ell} represents the transformation matrix for the
%%% nonlinearities of layer \ell to the base. To see what is base and what is transformation matrix I refer you to read the report on the Lip method shared in overleaf.
%%% The following lines creates these transformation matrices based on
%%% linear connections  between elements of base.
E_Q{1} = sparse(zeros(2*X(1), LMI_size));
index = system_dim;
E_Q{1}(X(1)+1:end  , index+1:index+X(1) ) =  sparse(eye(X(1)));
index = 0;
cntrb = W_n_n{1};
E_Q{1}(1:X(1), index+1:index+system_dim ) =  cntrb;


E_Q{2} = sparse(zeros(2*X(2), LMI_size));
index = system_dim+sum(X(1:1));
E_Q{2}(X(2)+1:end  , index+1:index+X(2) ) =  eye(X(2));
index = system_dim;
cntrb = W_n_n{2};
E_Q{2}(1:X(2), index+1:index+X(1) ) =  cntrb;
index = 0;
middle = 1;
cntrb = W_L_n{2}*middle*W_n_L{1};
E_Q{2}(1:X(2), index+1:index+system_dim ) =  cntrb;

for i = 3:len
    E_Q{i} = sparse(zeros(2*X(i), LMI_size));
    index = system_dim+sum(X(1:i-1));
    E_Q{i}(X(i)+1:end  , index+1:index+X(i) ) =  eye(X(i));
    index = system_dim+sum(X(1:i-1-1));
    cntrb = W_n_n{i};
    E_Q{i}(1:X(i), index+1:index+X(i-1) ) =  cntrb;
    index = system_dim+sum(X(1:i-2-1));
    middle = 1;
    cntrb = W_L_n{i}*middle*W_n_L{i-1};
    E_Q{i}(1:X(i), index+1:index+X(i-2) ) =  cntrb;
    for j = 3:i-1
        index = system_dim+sum(X(1:i-j-1));
        middle = middle*W_L_L{i-j+2};
        cntrb = W_L_n{i}*middle*W_n_L{i-j+1};
        E_Q{i}(1:X(i), index+1:index+X(i-j) ) =  cntrb;
    end
    index = 0;
    middle = middle*W_L_L{2};
    cntrb = W_L_n{i}*middle*W_n_L{1};
    E_Q{i}(1:X(i), index+1:index+system_dim ) =  cntrb;
end


%%%%  E_J represents the transformation matrix for Quadratic constraint
%%%%  which is written over the Lipstchitz inquality.
E_J = sparse(zeros(system_dim+1,LMI_size));
E_J(1:system_dim  ,  1:system_dim) = eye(system_dim);
E_J(system_dim+1,end-3:end) =  Net.weights{end};

%%%%%%%%%%%



constraints = [];   


%%% Here we capture three approaches for Lipstchitz constant computation:
%%%  1- LipSDP-netwrok  : In this technique we assign only one decision variable to all the Trapezium-FFNN. This technique sacrifises the  accuracy for a fast computation. But sometimes can also return  infeasibility as a singlr decision variable sometimes is not able to provide feasibly
%%%  2- LipsSDP-layer   : In this technique we assign one decision variable to every single layer in Trapezium-FFNN. This is slower than Lip-SDP network but result in more accurate Lip constant since contains more decision variables.
%%%  3- LipSDP-neurn    : This one is slowest but most accurate technique,it assigns a decision variable to every singlr nonlinear activation function. 

if strcmp(LipSDP_type, 'LipSDP-network')
    lambda = sdpvar(1,1);
    constraints = [constraints, lambda>=0];
end
for i = 1:len
    
    if strcmp(LipSDP_type, 'LipSDP-neuron')
        lambda{i} = sdpvar(X(i),1);
        for jj = 1:X(i)
            if ~strcmp(Net.layers{i}(jj+nl_index(i)-1), 'poslin' )  &&  ~strcmp(Net.layers{i}(jj+nl_index(i)-1), 'purelin' )
                
                constraints = [constraints, lambda{i}(jj,1)>=0];
                
            elseif strcmp(Net.layers{i}(jj+nl_index(i)-1), 'purelin' )
                
                error('Purelin neurons are included in the verification process')
                
            elseif strcmp(Net.layers{i}(jj+nl_index(i)-1), 'poslin' )
                
                if Lb{i}(jj+nl_index(i)-1 , 1)<0 && Ub{i}(jj+nl_index(i)-1 , 1)>0   %%% for a relu activation function if the function is always active/inactive then there is no reason to keep the decision variable positive, see the report. 
                    constraints = [constraints, lambda{i}(jj,1)>=0];
                end
                
            end
        end
        Q{i} = [-2*diag(alpha_params{i}.*beta_params{i})*diag(lambda{i})  diag((alpha_params{i}+beta_params{i}))*diag(lambda{i})   ;...
                 diag((alpha_params{i}+beta_params{i}))*diag(lambda{i})                      -2*diag(lambda{i})                    ];   %%%This is the quadratic constraint we use for Lip-SDP method, take a look at report. 
    
    elseif strcmp(LipSDP_type, 'LipSDP-layer')
        lambda{i} = sdpvar(1,1);
        constraints = [constraints, lambda{i}>=0];
        
        Q{i} = [-2*diag(alpha_params{i}.*beta_params{i})*lambda{i}  diag((alpha_params{i}+beta_params{i}))*lambda{i}  ;...
                 diag((alpha_params{i}+beta_params{i}))*lambda{i}                -2*lambda{i}*eye(X(i))               ];
    elseif strcmp(LipSDP_type, 'LipSDP-network')
   
        Q{i} = [-2*diag(alpha_params{i}.*beta_params{i})*lambda     diag((alpha_params{i}+beta_params{i}))*lambda     ;...
                 diag((alpha_params{i}+beta_params{i}))*lambda                      -2*lambda*eye(X(i))               ];
    end        
end 

rho = sdpvar(1,1);
constraints = [constraints, rho>=0];   %%% This is the main variable of optimization, this is the square of the local Lipstchitz constant for J_STL

%%% start of generation of our LMI
%%% This the part we add linear informations to the original quadratic constraint to make it feasble. see the report and study the computation of M_sdp.
M_sdp =  - E_J.' *[    rho*eye(system_dim)  zeros(system_dim,1)    ;...
                       zeros(1,system_dim)          -1             ]* E_J;

for t = 1:len
        M_sdp = M_sdp + E_Q{t}'*Q{t}*E_Q{t};
end

constraints = [constraints,M_sdp<=0];

%%% End of generation of our LMI

obj =  rho;

% yalmip_options = sdpsettings('solver','mosek','verbose',1, 'mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP',1);
yalmip_options = sdpsettings('solver','mosek','verbose',1);
out = optimize(constraints,obj,yalmip_options);
clc
time =  out.solvertime;
status  =  out.problem;
Lip    = double(sqrt(rho));

end