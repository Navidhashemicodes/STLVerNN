function [Lb, Ub, alpha_params, beta_params]= Trapezius_bounds_General(Center, epsilon, net, border,  method_Trajectory, method_STL)

%%% This function gets the set of intial condition (Center, epsilon) and
%%% introduce it to the trapezium-FFNN and then computes the pre-activation
%%% bounds. border is the border between the Trajectory and STL-FFNN in the
%%% trapezium-FFNN. Her we provide three techniques to compute
%%% pre-activation bounds on the Trajectory  
%%%  1 - method_Trajectory = Linear-CROWN in this method we utilize the linear technique proposed in CROWN library to compute the pre-activation bounds.
%%%  2 - method_Trajectory = Quadratic-CROWN in this method the CROWN library utilize lower and upper quadratic bounds on the activation function  and presents a convex quadratic programming through a couple of simplifications to compute the pre-activation bounds.
%%%  3 - method_Trajectory = Quadratic-Gurobi in this technique we again use the quadratic bounds on the activation functions to compute pre-activation bounds but we do not apply simplifications and present the problem as non-convex quadratic programming. 
     %%% Although this optimization is not convex but the gurobi guarantees they can provide us global minimum. This nonconvex quadratic programming is utilized to compute the pre-activation bounds.
%%% pre-activation bounds on the STL-FFNN : We know STL-FFNN includes only relu as its nonlinear activation functions. This implies we can utilize exact-star or approx-star technique to comute the preactivation bounds.
%%% 1 - method_STL= exact-star   ----> computes the pre-activation bounds with exact star technique.
%%% 2 - method_STL=approx-star ------> computes the pre-activation bounds with approx-star technique.

Net=Linear_Nonlinear(net, border);   %%% We utilize this function to make it certain the location of linear and nonlinear activations are completely separated in the Trapezium-FFNN.

switch method_Trajectory
    
    
    case 'Linear-CROWN'
        l0 = Center-epsilon;
        u0 = Center+epsilon;
        W_pos = 0.5*(Net.weights{1}+abs(Net.weights{1}));
        W_neg = 0.5*(Net.weights{1}-abs(Net.weights{1}));
        Lb{1} = W_pos*l0+W_neg*u0+Net.biases{1};
        Ub{1} = W_pos*u0+W_neg*l0+Net.biases{1};
        [param_l_linear{1}, param_u_linear{1} ] = linear_lines( Lb{1}, Ub{1}, net.layers{1} );
        [alpha_params{1},beta_params{1}] =   piecewise_linear_slopes(Lb{1}, Ub{1}, Net.layers{1});
        for i = 1:border-1
            disp(i)
            [Lb{i+1}, Ub{i+1}] = Linear_CROWN(net, i+1 , l0, u0, param_l_linear, param_u_linear);
            [param_l_linear{i+1}, param_u_linear{i+1} ] = linear_lines( Lb{i+1}, Ub{i+1}, net.layers{i+1} );
            [alpha_params{i+1},beta_params{i+1}] =   piecewise_linear_slopes(Lb{i+1}, Ub{i+1}, Net.layers{i+1});
        end
        
        STL_init = Star();
        STL_init.V =  [0.5*(Lb{end}+Ub{end}) eye(length(net.layers{border}))];
        STL_init.C =  zeros(1,length(net.layers{border}));
        STL_init.d =  0;
        STL_init.predicate_lb = -0.5*(Ub{end}-Lb{end});
        STL_init.predicate_ub =  0.5*(Ub{end}-Lb{end});
        STL_init.dim =  length(net.layers{border});
        STL.nVar =  length(net.layers{border});
        
        
        
        [Lb_STL, Ub_STL]  =  star_set(STL_init, border, length(Net.layers), Net, method_STL);  %%%the pre-activation bounds on the Trajectory is computed and this function computes pre-activation bounds on the STL-FFNN. 
        Lb(border+1:length(Net.layers)) = Lb_STL;
        Ub(border+1:length(Net.layers)) = Ub_STL;
        for i = border+1:length(Net.layers)
            [alpha_params{i} ,beta_params{i}]  =    piecewise_linear_slopes(Lb{i}, Ub{i}, Net.layers{i});
        end
        
        
    case 'Quadratic-CROWN'
        
        
        l0 = Center-epsilon;
        u0 = Center+epsilon;
        figure
        for i = 1:border-1
            L = net.layers{i};
            if strcmp(char(L(1)), 'purelin')
                nln_firstloc = 0;
                dd = true;
            else
                nln_firstloc = 1;
                dd = false;
            end
            
            while dd
                nln_firstloc = nln_firstloc+1;
                if ~strcmp(char(L(nln_firstloc)), 'purelin')
                    dd = false;
                end
            end
            %%%%% nln_first loc represents the first index for nonlinear
            %%%%% activation function
            nln_border{i} = 0;
            dd = true;
            while dd && (nln_firstloc>=2) && i<(border-1)   %%% This while loop searchs for the independent linear layers and removes them from pre-activation bound computation. They are linear and their pre-activation bounds are obvious.
                SS1 =  zeros(   nln_firstloc-1  ,  size(net.weights{i},2)  );   SS1(1:nln_firstloc-1, 1:nln_firstloc-1) = eye(nln_firstloc-1); 
                SS2 =  zeros(   size(net.weights{i},1)  ,  nln_firstloc-1  );   SS2(1:nln_firstloc-1, 1:nln_firstloc-1) = eye(nln_firstloc-1);
                if nln_firstloc>1 && logical(prod(prod(net.biases{i}(1:nln_firstloc-1)==zeros(nln_firstloc-1,1)))) && logical(prod(prod(net.weights{i}(1:nln_firstloc-1,:)==SS1)))...
                                  && logical(prod(prod(net.weights{i}(:, 1:nln_firstloc-1)==SS2)))
                    SS1 =  zeros(   nln_firstloc-1  ,  size(net.weights{i+1},2)  );   SS1(1:nln_firstloc-1, 1:nln_firstloc-1) = eye(nln_firstloc-1); 
                    SS2 =  zeros(   size(net.weights{i+1},1)  ,  nln_firstloc-1  );   SS2(1:nln_firstloc-1, 1:nln_firstloc-1) = eye(nln_firstloc-1);
                    if logical(prod(prod(net.biases{i+1}(1:nln_firstloc-1)==zeros(nln_firstloc-1,1)))) && logical(prod(prod(net.weights{i+1}(1:nln_firstloc-1,:)==SS1)))...
                                                                                                       && logical(prod(prod(net.weights{i+1}(:, 1:nln_firstloc-1)==SS2)))
                        nln_border{i} = nln_firstloc-1;
                        dd = false;
                    else
                        nln_firstloc = nln_firstloc-1;
                    end
                else
                    nln_firstloc = nln_firstloc-1;
                end
            end
            plot(i,nln_border{i},'*')
            hold on
            %%%% We are pretty sure that from 1 upto nln_border the network
            %%%% is identity and independent
            if i==1
                index = nln_border{1};
                W1 = Net.weights{1}(index+1:end, index+1:end);
                b1 = Net.biases{1}(index+1:end, 1);
                W2 = Net.weights{2}(index+1:end, index+1:end);
                b2 = Net.biases{2}(index+1:end, 1);
                W_pos = 0.5*(Net.weights{1}+abs(Net.weights{1}));
                W_neg = 0.5*(Net.weights{1}-abs(Net.weights{1}));
                Lb{1} = W_pos*l0+W_neg*u0+Net.biases{1};
                Ub{1} = W_pos*u0+W_neg*l0+Net.biases{1};
                [param_l,param_u]   =    quad_low_and_up(Lb{1}(index+1:end), Ub{1}(index+1:end), Net.layers{1}(index+1:end ,1));  %%% returns the best upper and lower quadratic bounds on a single activation function
                [param_l_linear, param_u_linear ] = linear_lines( Lb{1}(index+1:end), Ub{1}(index+1:end), Net.layers{1}(index+1:end ,1) );  %%% whenever the quadratic bound provides non-convexity the CROWN library
                                                                                                                                            %%% replaces it with linear bound. Here the linear lines is in fact providing the upper and lower linear bounds on a single activation function
                [alpha_params{1},beta_params{1}] =   piecewise_linear_slopes(Lb{1}, Ub{1}, Net.layers{1});   %%% Given the pre-activation bounds on the layer, this function generates the bounds over the slope of activation function.
                l_init = l0;
                u_init = u0;
                [LLb{2}, UUb{2}] = Convex_quad_prebound(  W1, b1, W2, b2, l_init, u_init, param_l, param_u, param_u_linear, param_l_linear  );  %%% This function generates the convex quadratic programming and optimizes the pre-activation bounds.
                [alpha_params{2} ,beta_params{2}] =   piecewise_linear_slopes(LLb{2}, UUb{2}, Net.layers{2});
                Lb{2} = [Lb{1}(1:index,1);LLb{2}];
                Ub{2} = [Ub{1}(1:index,1);UUb{2}];
            else
                index = nln_border{i};
                W1 = Net.weights{i}(index+1:end, index+1:end);
                b1 = Net.biases{i}(index+1:end, 1);
                W2 = Net.weights{i+1}(index+1:end, index+1:end);
                b2 = Net.biases{i+1}(index+1:end, 1);
                l_init = [];
                u_init = [];
                L = Net.layers{i-1};
                for ii = index+1:length(L)
                    func = char(L(ii));
                    funci = str2func(func);
                    l_init = [l_init; funci(Lb{i-1}(ii,1)) ];   u_init = [u_init; funci(Ub{i-1}(ii,1)) ];
                end
                [param_l,param_u]   =    quad_low_and_up(Lb{i}(index+1:end), Ub{i}(index+1:end), Net.layers{i}(index+1:end ,1));
                [param_l_linear, param_u_linear ] = linear_lines( Lb{i}(index+1:end), Ub{i}(index+1:end), Net.layers{i}(index+1:end ,1) );
                [LLb{i+1}, UUb{i+1}] = Convex_quad_prebound(  W1, b1, W2, b2, l_init, u_init, param_l, param_u, param_u_linear, param_l_linear  );
                [alpha_params{i+1} ,beta_params{i+1}] =   piecewise_linear_slopes(LLb{i+1}, UUb{i+1}, Net.layers{i+1});
                Lb{i+1} = [Lb{i}(1:index,1);LLb{i+1}];
                Ub{i+1} = [Ub{i}(1:index,1);UUb{i+1}];
            end
        end
        
        
        STL_init = Star();
        STL_init.V =  [0.5*(Lb{end}+Ub{end}) eye(length(net.layers{border}))];
        STL_init.C =  zeros(1,length(net.layers{border}));
        STL_init.d =  0;
        STL_init.predicate_lb = -0.5*(Ub{end}-Lb{end});
        STL_init.predicate_ub =  0.5*(Ub{end}-Lb{end});
        STL_init.dim =  length(net.layers{border});
        STL.nVar =  length(net.layers{border});
        
        
        
        [Lb_STL, Ub_STL] = star_set(STL_init, border, length(Net.layers), Net, method_STL);  %%%the pre-activation bounds on the Trajectory is computed and this function computes pre-activation bounds on the STL-FFNN.
        Lb(border+1:length(Net.layers)) = Lb_STL;
        Ub(border+1:length(Net.layers)) = Ub_STL;
        for i = border+1:length(Net.layers)
            [alpha_params{i} ,beta_params{i}] =   piecewise_linear_slopes(Lb{i}, Ub{i}, Net.layers{i});
        end
        
        
        
        
    case 'Quadratic-Gurobi'
        
        
        l0 = Center-epsilon;
        u0 = Center+epsilon;
        for i = 1:border-1
            L = net.layers{i};
            if strcmp(char(L(1)), 'purelin')
                nln_firstloc = 0;
                dd = true;
            else
                nln_firstloc = 1;
                dd = false;
            end
            
            while dd
                nln_firstloc = nln_firstloc+1;
                if ~strcmp(char(L(nln_firstloc)), 'purelin')
                    dd = false;
                end
            end
            %%%%% nln_first loc represents the first index for nonlinear
            %%%%% activation function
            nln_border{i} = 0;
            dd = true;
            while dd && (nln_firstloc>=2) && i<(border-1)
                SS1 =  zeros(   nln_firstloc-1  ,  size(net.weights{i},2)  );   SS1(1:nln_firstloc-1, 1:nln_firstloc-1) = eye(nln_firstloc-1); 
                SS2 =  zeros(   size(net.weights{i},1)  ,  nln_firstloc-1  );   SS2(1:nln_firstloc-1, 1:nln_firstloc-1) = eye(nln_firstloc-1);
                if nln_firstloc>1 && logical(prod(prod(net.biases{i}(1:nln_firstloc-1)==zeros(nln_firstloc-1,1)))) && logical(prod(prod(net.weights{i}(1:nln_firstloc-1,:)==SS1)))...
                                  && logical(prod(prod(net.weights{i}(:, 1:nln_firstloc-1)==SS2)))
                    SS1 =  zeros(   nln_firstloc-1  ,  size(net.weights{i+1},2)  );   SS1(1:nln_firstloc-1, 1:nln_firstloc-1) = eye(nln_firstloc-1); 
                    SS2 =  zeros(   size(net.weights{i+1},1)  ,  nln_firstloc-1  );   SS2(1:nln_firstloc-1, 1:nln_firstloc-1) = eye(nln_firstloc-1);
                    if logical(prod(prod(net.biases{i+1}(1:nln_firstloc-1) ==zeros(nln_firstloc-1,1)))) && logical(prod(prod(net.weights{i+1}(1:nln_firstloc-1,:)==SS1)))...
                                                                                                       && logical(prod(prod(net.weights{i+1}(:, 1:nln_firstloc-1)==SS2)))
                        nln_border{i} = nln_firstloc-1;
                        dd = false;
                    else
                        nln_firstloc = nln_firstloc-1;
                    end
                else
                    nln_firstloc = nln_firstloc-1;
                end
            end
            
            %%%% We are pretty sure that from1 upto nln_border the network
            %%%% is identity and independent
            if i==1
                index = nln_border{1};
                W1 = Net.weights{1}(index+1:end, index+1:end);
                b1 = Net.biases{1}(index+1:end, 1);
                W2 = Net.weights{2}(index+1:end, index+1:end);
                b2 = Net.biases{2}(index+1:end, 1);
                W_pos = 0.5*(Net.weights{1}+abs(Net.weights{1}));
                W_neg = 0.5*(Net.weights{1}-abs(Net.weights{1}));
                Lb{1} = W_pos*l0+W_neg*u0+Net.biases{1};
                Ub{1} = W_pos*u0+W_neg*l0+Net.biases{1};
                [alpha_params{1} ,beta_params{1}] =   piecewise_linear_slopes(Lb{1}, Ub{1}, Net.layers{1});
                [param_l{1},param_u{1}]   =    quad_low_and_up(Lb{1}(index+1:end), Ub{1}(index+1:end), Net.layers{1}(index+1:end ,1));
                l_init = l0;
                u_init = u0;
                [LLb{2}, UUb{2}] = Gurobi_quad_prebound(  W1, b1, W2, b2, l_init, u_init, param_l{1}, param_u{1} ); 
                [alpha_params{2} ,beta_params{2}] =   piecewise_linear_slopes(LLb{2}, UUb{2}, Net.layers{2});
                Lb{2} = [Lb{1}(1:index,1);LLb{2}];
                Ub{2} = [Ub{1}(1:index,1);UUb{2}];
            else
                index = nln_border{i};
                W1 = Net.weights{i}(index+1:end, index+1:end);
                b1 = Net.biases{i}(index+1:end, 1);
                W2 = Net.weights{i+1}(index+1:end, index+1:end);
                b2 = Net.biases{i+1}(index+1:end, 1);
                l_init = [];
                u_init = [];
                L = Net.layers{i-1};
                for ii = index+1:length(L)
                    func = char(L(ii));
                    funci = str2func(func);
                    l_init = [l_init; funci(Lb{i-1}(ii,1)) ]  ; u_init = [u_init; funci(Ub{i-1}(ii,1)) ]  ;
                end
                [param_l{i},param_u{i}]   =    quad_low_and_up(Lb{i}(index+1:end), Ub{i}(index+1:end), Net.layers{i}(index+1:end ,1));
                [LLb{i+1}, UUb{i+1}] = Gurobi_quad_prebound(  W1, b1, W2, b2, l_init, u_init, param_l{i}, param_u{i} );
                [alpha_params{i+1} ,beta_params{i+1}] =   piecewise_linear_slopes(LLb{i+1}, UUb{i+1}, Net.layers{i+1});
                Lb{i+1} = [Lb{i}(1:index,1);LLb{i+1}];
                Ub{i+1} = [Ub{i}(1:index,1);UUb{i+1}];
            end
        end
        
        STL_init = Star();
        STL_init.V =  [0.5*(Lb{end}+Ub{end}) eye(length(net.layers{border}))];
        STL_init.C =  zeros(1,length(net.layers{border}));
        STL_init.d =  0;
        STL_init.predicate_lb = -0.5*(Ub{end}-Lb{end});
        STL_init.predicate_ub =  0.5*(Ub{end}-Lb{end});
        STL_init.dim =  length(net.layers{border});
        STL.nVar =  length(net.layers{border});
        
        
        
        [Lb_STL, Ub_STL] = star_set(STL_init, border, length(Net.layers), Net, method_STL);
        Lb(border+1:length(Net.layers)) = Lb_STL;
        Ub(border+1:length(Net.layers)) = Ub_STL;
        for i = border+1:length(Net.layers)
            [alpha_params{i} ,beta_params{i}] =   piecewise_linear_slopes(Lb{i}, Ub{i}, Net.layers{i});
        end
        
        
end

end




function [LLb, UUb] = star_set(In, start, the_end, Net, analysis_type)   %%% receives the input of STL-FFNN as a star set and returns the pre-activation bounds on every single layer for STL-FFNN
in{1} = In;
for jj = 1:the_end-start+1
    i = jj+start-1;

    disp(['exact reachability analysis for Layer ' num2str(jj) ' of ' num2str(the_end-start+1) ' in STL ReLU network.'])
    lenS = length(in{jj});
    lenAC = 0;
    for k = 1:lenS
        if jj ==1
            in_layer{jj}(k) = in{jj}(k);
        else
            in_layer{jj}(k) = affineMap(in{jj}(k), Net.weights{i}, Net.biases{i});
        end
        In_layer = in_layer{jj}(k);
        if strcmp(Net.layers{i}(1), 'purelin')
            index = 1;
            dd = true;
            while dd
                index = index+1;
                if ~strcmp(Net.layers{i}(index), 'purelin')
                    dd = false;
                end
            end
        else
            index = 1;
        end
        In_layer.V = In_layer.V(index:end, :);
        if size(In_layer.state_lb , 1 )~= 0
            In_layer.state_lb = In_layer.state_lb(index:end, 1);
            In_layer.state_ub = In_layer.state_ub(index:end, 1);
        end
        In_layer.dim =  size(In_layer.V , 1 );
        
        
        Layers = LayerS(eye(In_layer.dim), zeros(In_layer.dim,1) , 'poslin');
        Ln     = LayerS(eye(In_layer.dim), zeros(In_layer.dim,1) , 'purelin');
        Layers = [Layers Ln];
        
        F = FFNNS(Layers);
        
        [Out_layer, ~] = F.reach(In_layer, analysis_type);
        
        lenAC_now = length(Out_layer);
        for j = 1:lenAC_now
            if strcmp(analysis_type, 'exact-star')
                Out_layer(j).V = [in_layer{jj}(k).V(1:index-1,:)  ;  Out_layer(j).V];   %%% you should check if the predicate updating process affects the V or not, I guess No, but you need to check
            elseif strcmp(analysis_type, 'approx-star')
                d_size = size(Out_layer(j).V,2)-size(in_layer{jj}(k).V(1:index-1,:), 2);
                Out_layer(j).V = [[in_layer{jj}(k).V(1:index-1,:) zeros(index-1, d_size) ] ;  Out_layer(j).V];
            end
            
            if size(in_layer{jj}(k).state_lb , 1 )~=0
                if size(Out_layer(j).state_lb , 1 )~=0
                    Out_layer(j).state_lb = [in_layer{jj}(k).state_lb(1:index-1,1) ; Out_layer(j).state_lb];
                    Out_layer(j).state_ub = [in_layer{jj}(k).state_ub(1:index-1,1) ; Out_layer(j).state_ub];
                end
            else
                if size(Out_layer(j).state_lb , 1 )~=0
                    SSS =  getBox(in_layer{jj}(k));
                    Out_layer(j).state_lb = [SSS.lb(1:index-1,1) ; Out_layer(j).state_lb];
                    Out_layer(j).state_ub = [SSS.ub(1:index-1,1) ; Out_layer(j).state_ub];
                end
            end

            Out_layer(j).dim = Out_layer(j).dim + (index-1);
            in{jj+1}(lenAC+j) = Out_layer(j);
        end
        lenAC     = lenAC+ lenAC_now;
    end 
end

disp('pre-activation computation is in process ...')
%%%%%%  returns the pre-activation bounds over the layers  please remove
%%%%%%  the pre-activation computations for purelins later.
for jj = 1:the_end-start+1
    i = jj+start-1;
    clear the_lb the_ub
    for j = 1:length(in_layer{jj})
        if size(in_layer{jj}(j).state_lb , 1 ) ==0
            SSS = getBox(in_layer{jj}(j));
            the_lb(:,j) = SSS.lb;
            the_ub(:,j) = SSS.ub;
        else
            the_lb(:,j) = in_layer{jj}(j).state_lb;
            the_ub(:,j) = in_layer{jj}(j).state_ub;
        end
    end
    if size(the_lb,2)>1
        Lb{jj} = min(the_lb')';
        Ub{jj} = max(the_ub')';
    else
        Lb{jj} = the_lb;
        Ub{jj} = the_ub;
    end
end
LLb = Lb(2:end);
UUb = Ub(2:end);
end




function [l,u] = Gurobi_quad_prebound( W_pre, b_pre, W, b, l_init, u_init, param_u, param_l ) 

%%% This function receives
%%%    1- the upper and lower quadratic functions (param_u, param_l) on each activation of the layer, 
%%%    2- it also receives the Weight and bias vector before the layer (W_pre, b_pre)
%%%    3- and it utilizes the pre-activation bounds on the previous layer (l_init, u_init)
%%%    4- and the weights and bias vector after the layer (W, b),
%%% to characterize the a qudratic non-convex programming and solves it to compute the pre-activation bounds in the current layer. This non-convex optimization will be handled by Gurobi solvers (bilinear programming).

len1 = size(b,1);
len2 = size(b_pre,1);
len3 = size(W_pre,2);
clear model
for j = 1:len1
    q_U{j} = zeros(len2,1);
    Lambda_U{j} = zeros(1,len2);
    Sigma_U{j} = zeros(len2,1);
    for i = 1:len2
        if W(j,i)>0
            q_U{j}(i,1) = W(j,i)*param_u(1,i);
            Lambda_U{j}(1,i) = W(j,i)*param_u(2,i);
            Sigma_U{j}(i,1) = param_u(3,i);
        else
            q_U{j}(i,1) = W(j,i)*param_l(1,i);
            Lambda_U{j}(1,i) = W(j,i)*param_l(2,i);
            Sigma_U{j}(i,1) = param_l(3,i);
        end
    end
    
    cd C:\gurobi951\win64\matlab
    model.A =  sparse(len3,len3);
    model.Q  =  sparse( W_pre.' * diag(q_U{j}) * W_pre );
    model.obj  =  2 * b_pre.' * diag(q_U{j}) * W_pre  + Lambda_U{j} * W_pre;
    model.lb =  l_init;
    model.ub =   u_init;
    model.modelsense =  'Max';
    params.NonConvex  =  2;
    gurobi_write(model, 'qp.lp',params);
    results = gurobi(model,params);
    
    u(j,1) =  results.objval  +   b_pre.' * diag(q_U{j}) * b_pre  +  Lambda_U{j} * b_pre  +  W(j,:) * Sigma_U{j} + b(j,1);
    cd C:\Users\Navid\Documents\MATLAB\Others\Personal_Trapezius_toolbox
end



clear model
for j = 1:len1
    q_L{j} = zeros(len2,1);
    Lambda_L{j} = zeros(1,len2);
    Sigma_L{j} = zeros(len2,1);
    for i = 1:len2
        if W(j,i)>0
            q_L{j}(i,1) = W(j,i)*param_l(1,i);
            Lambda_L{j}(1,i) = W(j,i)*param_l(2,i);
            Sigma_L{j}(i,1) = param_l(3,i);
        else
            q_L{j}(i,1) = W(j,i)*param_u(1,i);
            Lambda_L{j}(1,i) = W(j,i)*param_u(2,i);
            Sigma_L{j}(i,1) = param_u(3,i);
        end
    end
    
    cd C:\gurobi951\win64\matlab
    model.A =  sparse(len3,len3);
    model.Q  =  sparse( W_pre.' * diag(q_L{j}) * W_pre );
    model.obj  =  2 * b_pre.' * diag(q_L{j}) * W_pre  + Lambda_L{j} * W_pre;
    model.lb =  l_init;
    model.ub =   u_init;
    model.params.NonConvex = 2;
    model.modelsense =  'Min';
    params.NonConvex = 2;
    gurobi_write(model, 'qp.lp',params);
    results = gurobi(model,params);
    
    l(j,1) =  results.objval  +   b_pre.' * diag(q_L{j}) * b_pre  +  Lambda_L{j} * b_pre  +  W(j,:) * Sigma_L{j} + b(j,1);
    cd C:\Users\Navid\Documents\MATLAB\Others\Personal_Trapezius_toolbox
end
    
end






function [l,u] = Convex_quad_prebound( W_pre, b_pre, W, b, l_init, u_init, param_u, param_l, param_u_linear, param_l_linear  )

%%% This function receives
%%%    1- the upper and lower quadratic functions (param_u, param_l) on each activation of the layer,
%%%    2- the upper and lower linear functions (param_u_linear, param_l_linear) on each activation of the layer,
%%%    3- it also receives the Weight and bias vector before the layer (W_pre, b_pre)
%%%    4- and it utilizes the pre-activation bounds on the previous layer (l_init, u_init)
%%%    5- and the weights and bias vector after the layer (W, b),
%%% to characterize the a qudratic convex programming and solves it to compute the pre-activation bounds in the current layer. 

len1 = size(b,1);
len2 = size(b_pre,1);
len3 = size(W_pre,2);
clear model
for j = 1:len1
    q_U{j} = zeros(len2,1);
    Lambda_U{j} = zeros(1,len2);
    Sigma_U{j} = zeros(len2,1);
    for i = 1:len2
        if W(j,i)>0
            if param_u(1,i)<0
                q_U{j}(i,1) = W(j,i)*param_u(1,i);
                Lambda_U{j}(1,i) = W(j,i)*param_u(2,i);
                Sigma_U{j}(i,1) = param_u(3,i);
            else
                q_U{j}(i,1) = 0;
                Lambda_U{j}(1,i) = W(j,i)*param_u_linear(1,i);
                Sigma_U{j}(i,1) = param_u_linear(2,i);
            end
        else
            if param_l(1,i)>0
                q_U{j}(i,1) = W(j,i)*param_l(1,i);
                Lambda_U{j}(1,i) = W(j,i)*param_l(2,i);
                Sigma_U{j}(i,1) = param_l(3,i);
            else
                q_U{j}(i,1) = 0;
                Lambda_U{j}(1,i) = W(j,i)*param_l_linear(1,i);
                Sigma_U{j}(i,1) = param_l_linear(2,i);
            end
        end
    end
    
    cd C:\gurobi951\win64\matlab
    model.A =  sparse(len3,len3);
    model.Q  =  sparse( W_pre.' * diag(q_U{j}) * W_pre );
    model.obj = 2 * b_pre.' * diag(q_U{j}) * W_pre  + Lambda_U{j} * W_pre;
    model.lb =  l_init;
    model.ub =   u_init;
    model.modelsense =  'Max';
    gurobi_write(model, 'qp.lp');
    results  =  gurobi(model);
    
    u(j,1) =  results.objval  +   b_pre.' * diag(q_U{j}) * b_pre  +  Lambda_U{j} * b_pre  +  W(j,:) * Sigma_U{j} + b(j,1);
    cd C:\Users\Navid\Documents\MATLAB\Others\Personal_Trapezius_toolbox
end



clear model
for j = 1:len1
    q_L{j} = zeros(len2,1);
    Lambda_L{j} = zeros(1,len2);
    Sigma_L{j} = zeros(len2,1);
    for i = 1:len2
        if W(j,i)>0
            if param_l(1,i)>0
                q_L{j}(i,1) = W(j,i)*param_l(1,i);
                Lambda_L{j}(1,i) = W(j,i)*param_l(2,i);
                Sigma_L{j}(i,1) = param_l(3,i);
            else
                q_L{j}(i,1) = 0;
                Lambda_L{j}(1,i) = W(j,i)*param_l_linear(1,i);
                Sigma_L{j}(i,1) = param_l_linear(2,i);
            end
        else
            if param_u(1,i)<0
                q_L{j}(i,1) = W(j,i)*param_u(1,i);
                Lambda_L{j}(1,i) = W(j,i)*param_u(2,i);
                Sigma_L{j}(i,1) = param_u(3,i);
            else
                q_L{j}(i,1) = 0;
                Lambda_L{j}(1,i) = W(j,i)*param_u_linear(1,i);
                Sigma_L{j}(i,1) = param_u_linear(2,i);
            end
        end
    end
    
    cd C:\gurobi951\win64\matlab
    model.A =  sparse(len3,len3);
    model.Q = sparse( W_pre.' * diag(q_L{j}) * W_pre );
    model.obj = 2 * b_pre.' * diag(q_L{j}) * W_pre  + Lambda_L{j} * W_pre;
    model.lb =  l_init;
    model.ub =   u_init;
    model.modelsense =  'Min';
    gurobi_write(model, 'qp.lp');
    results = gurobi(model);
    
    l(j,1) =  results.objval  +   b_pre.' * diag(q_L{j}) * b_pre  +  Lambda_L{j} * b_pre  +  W(j,:) * Sigma_L{j} + b(j,1);
    cd C:\Users\Navid\Documents\MATLAB\Others\Personal_Trapezius_toolbox
end
    
end






function [param_l,param_u]   =    quad_low_and_up(l,u, funcs)

%%% consider a layer od activation functions (funcs).
%%% This function receives the pre-activan bounds on this layer (l,u), and
%%% charactrizes the two quadratic functions for each signle activation.
%%% The former is always an upper bound and the latter is always lower
%%% bound.


NN = 1000;
for i = 1:length(l)
    if strcmp(funcs(i), 'poslin')
        
        if l(i)>0
            param_u(:,i) = [0;1;0];
            param_l(:,i) = [0;1;0];
        elseif u(i)<=0
            param_u(:,i) = [0;0;0];
            param_l(:,i) = [0;0;0];
        else
            if u(i)-l(i)>300
                param_u(:,i) = [0 ;u(i)/(u(i)-l(i)) ; -u(i)*l(i)/(u(i)-l(i)) ];
            else
                A = zeros(NN-2,3);
                b = zeros(NN-2,1);
                sample = linspace(l(i),u(i),NN);
                A_eq = [sample(1)^2, sample(1), 1; sample(end)^2, sample(end), 1];
                b_eq = [0;sample(end)];
                for j = 2:NN-1
                    A(j-1,:) = [sample(j)^2, sample(j), 1];
                    b(j-1,:) = f_x(sample(j));
                end
                cd C:\gurobi951\win64\examples\matlab
                f = mean(A);
                [param_u(:,i),~,~,~,~]  =  linprog(f,-A,-b,A_eq,b_eq,[],[],[],[]);
                cd   C:\Users\Navid\Documents\MATLAB\Others\Personal_Trapezius_toolbox
            end
            param_l(:,i) = [1/(u(i)-l(i))  ;  -l(i)/(u(i)-l(i)) ;  0];
        end
        
    elseif strcmp(funcs(i), 'purelin')
        param_u(:,i) = [0;1;0];
        param_l(:,i) = [0;1;0];
    else
        func = char(funcs(i));
        f_x = str2func(func);
        A = zeros(NN,3);
        b = zeros(NN,1);
        sample = linspace(l(i),u(i),NN);
        for j = 1:NN
            A(j,:) = [sample(j)^2, sample(j), 1];
            b(j,:) = f_x(sample(j));
        end
        cd C:\gurobi951\win64\examples\matlab
        f = mean(A);
        [param_u(:,i),~,~,~,~]  =  linprog(f,-A,-b,[],[],[],[],[],[]);
        f = -mean(A);
        [param_l(:,i),~,~,~,~] = linprog(f, A, b,[],[],[],[],[],[]);
        cd   C:\Users\Navid\Documents\MATLAB\Others\Personal_Trapezius_toolbox
    end
end
end



function   [param_l_linear, param_u_linear] = linear_lines(l_pre,u_pre, funcs)

%%% consider a layer of activation functions (funcs).
%%% This function receives the pre-activan bounds on this layer (l_pre,u_pre), and
%%% charactrizes two linear functions for each signle activation.
%%% The former is always an upper bound and the latter is always lower
%%% bound.

leng = length(l_pre);
alph_upp = [];
alph_low = [];
alphabet_upp = [];
alphabet_low = [];

for i = 1:leng
    
    if strcmp(funcs(i), 'poslin')
        if l_pre(i,1)>0  
            
            alph_upp =  1;      
            alphabet_upp =  0; 
            
            alph_low =  1 ;
            alphabet_low =  0  ;
            
        elseif u_pre(i,1)<0 
            
            alph_upp =  0 ;                 
            alphabet_upp =  0 ; 
            
            alph_low =  0 ;         
            alphabet_low =  0 ;  
            
        else 
            
            ss = (u_pre(i,1)-0)/(u_pre(i,1)-l_pre(i,1));
            alph_upp =  ss   ;
            alphabet_upp =  -ss*l_pre(i,1)  ;
            
            alph_low =  ss;
            alphabet_low =  0 ;
            
        end
    elseif strcmp(funcs(i), 'purelin')
        
            alph_upp = 1;      
            alphabet_upp = 0; 
            
            alph_low = 1 ;
            alphabet_low =  0  ;
        
    elseif strcmp(funcs(i), 'tanh')
        if l_pre(i,1)>0
            d = 0.5*(  l_pre(i,1) +   u_pre(i,1)   );
            alph_upp =  1-tanh(d)^2;
            alphabet_upp =   tanh(d)-(1-tanh(d)^2)*d;
            ss = (  tanh(u_pre(i,1))-tanh(l_pre(i,1))  ) / (    u_pre(i,1)-l_pre(i,1)   );
            alph_low =  ss ;
            alphabet_low =  tanh(l_pre(i,1))-ss*l_pre(i,1)  ;
            
        elseif u_pre(i,1)<0
            ss = (tanh(u_pre(i,1))-tanh(l_pre(i,1)))/(u_pre(i,1)-l_pre(i,1));
            alph_upp =  ss ;
            alphabet_upp =  tanh(l_pre(i,1))-ss*l_pre(i,1)   ;
            d = 0.5*(l_pre(i,1)+u_pre(i,1));
            alph_low =  1-tanh(d)^2 ;
            alphabet_low = tanh(d)-(1-tanh(d)^2)*d ;
        else
            x0 = 0.01;
            fun = @(x)d_find_tanh(x, l_pre(i,1));
%             option=optimset('MaxFunEvals',500000);
%             option=optimset(option,'MaxIter',500);
%             option=optimset(option,'disp','iter','LargeScale','off','TolFun',0.0001);
            [d,fval]=fmincon(fun,x0,[],[],[],[],l_pre(i,1),[],[]);
            disp(fval);
            if (d<0 || d>u_pre(i,1))
                ss = (tanh(u_pre(i,1))-tanh(l_pre(i,1)))/(u_pre(i,1)-l_pre(i,1));
                alph_upp = ss   ;
                alphabet_upp = tanh(l_pre(i,1))-ss*l_pre(i,1) ;
            else
                alph_upp = 1-tanh(d)^2;
                alphabet_upp = tanh(l_pre(i,1))-(1-tanh(d)^2)*l_pre(i,1) ;
            end
            
            x0 = 0.01;
            fun = @(x)d_find_tanh(x, u_pre(i,1));
%             option=optimset('MaxFunEvals',500000);
%             option=optimset(option,'MaxIter',500);
%             option=optimset(option,'disp','iter','LargeScale','off','TolFun',0.0001);
            [d,fval]=fmincon(fun,x0,[],[],[],[],[],u_pre(i,1),[]);
            disp(fval);
            if (d>0 || d<l_pre(i,1))
                
                ss = (tanh(u_pre(i,1))-tanh(l_pre(i,1)))/(u_pre(i,1)-l_pre(i,1));
                alph_low =  ss ;
                alphabet_low = tanh(u_pre(i,1))-ss*u_pre(i,1) ;
            else
                alph_low = 1-tanh(d)^2;
                alphabet_low =  tanh(u_pre(i,1))-(1-tanh(d)^2)*u_pre(i,1)  ;
            end
        end
    elseif strcmp(funcs(i), 'sigmoid')
        if l_pre(i,1)>0
            d = 0.5*(  l_pre(i,1) +   u_pre(i,1)   );
            alph_upp =  sigmoid(d)*(1-sigmoid(d));
            alphabet_upp =   sigmoid(d)-(sigmoid(d)*(1-sigmoid(d)))*d;
            ss = (  sigmoid(u_pre(i,1))-sigmoid(l_pre(i,1))  ) / (    u_pre(i,1)-l_pre(i,1)   );
            alph_low =  ss ;
            alphabet_low =  sigmoid(l_pre(i,1))-ss*l_pre(i,1)  ;
            
        elseif u_pre(i,1)<0
            ss = (sigmoid(u_pre(i,1))-sigmoid(l_pre(i,1)))/(u_pre(i,1)-l_pre(i,1));
            alph_upp =  ss ;
            alphabet_upp =  sigmoid(l_pre(i,1))-ss*l_pre(i,1)   ;
            d = 0.5*(l_pre(i,1)+u_pre(i,1));
            alph_low =  sigmoid(d)*(1-sigmoid(d)) ;
            alphabet_low = sigmoid(d)-sigmoid(d)*(1-sigmoid(d))*d ;
        else
            x0 = 0.01;
            fun = @(x)d_find_sigmoid(x, l_pre(i,1));
            option=optimset('MaxFunEvals',500000);
            option=optimset(option,'MaxIter',500);
            option=optimset(option,'disp','iter','LargeScale','off','TolFun',0.0001);
            [d,fval]=fsolve(fun,x0,option);
            disp(fval);
            clc
%             roots = fzero(fun,x0);
%             d = roots(1);
            if (d<0 || d>u_pre(i,1))
                ss = (sigmoid(u_pre(i,1))-sigmoid(l_pre(i,1)))/(u_pre(i,1)-l_pre(i,1));
                alph_upp = ss   ;
                alphabet_upp = sigmoid(l_pre(i,1))-ss*l_pre(i,1) ;
            else
                alph_upp = sigmoid(d)*(1-sigmoid(d));
                alphabet_upp = sigmoid(l_pre(i,1))-sigmoid(d)*(1-sigmoid(d))*l_pre(i,1) ;
            end
            
            x0 = 0.01;
            fun = @(x)d_find_sigmoid(x, u_pre(i,1));
            option=optimset('MaxFunEvals',500000);
            option=optimset(option,'MaxIter',500);
            option=optimset(option,'disp','iter','LargeScale','off','TolFun',0.0001);
            [d,fval]=fsolve(fun,x0,option);
            disp(fval);
            clc
%             roots = fzero(fun,x0);
%             d = roots(1);
            if (d>0 || d<l_pre(i,1))
                
                ss = (sigmoid(u_pre(i,1))-sigmoid(l_pre(i,1)))/(u_pre(i,1)-l_pre(i,1));
                alph_low =  ss ;
                alphabet_low = sigmoid(u_pre(i,1))-ss*u_pre(i,1) ;
            else
                alph_low = sigmoid(d)*(1-sigmoid(d));
                alphabet_low =  sigmoid(u_pre(i,1))-sigmoid(d)*(1-sigmoid(d))*u_pre(i,1)  ;
            end
        end
    else
        error('the mentioned activation function is not supported')
    end
    param_l_linear(:,i) = [alph_low;alphabet_low];
    param_u_linear(:,i) = [alph_upp;alphabet_upp];
end
end

function Eq  =  d_find_tanh(x, data)
    Eq = (((tanh(x)-tanh(data))/(x-data))-(1-tanh(x)^2))^2;
end

function Eq  =  d_find_sigmoid(x, data)
    Eq = (((sigmoid(x)-sigmoid(data))/(x-data))-sigmoid(x)*(1-sigmoid(x)))^2;
end

function   [Lbm, Ubm, param_l_linear, param_u_linear] = Linear_CROWN(net, m , l0, u0, param_l_linear, param_u_linear )

%%%% receives the linear bounds (param_l_linear, param_u_linear) for layers i=1,...,m-1 and the
%%%% initial state's bounds (l0, u0) and compute the linear bounds on the
%%%% input of layer "m" from neural network "net" (param_l_linear, param_u_linear). It also presents the
%%%% pre-activation bounds on the m-th layer (Lbm, Ubm).

n0 = size(net.weights{1},2);
for i = 1:m
    n{i} = length(net.layers{i});
end


I = eye(n{m});
for j = 1:n{m}
    
    Lambda{m}(j,:) = I(j,:);
    Omega{m}(j,:) = I(j,:);
    
    Delta{m}(:,j) = zeros(n{m-1},1);
    Theta{m}(:,j) = zeros(n{m-1},1);
    intersection_upp{j} =  Lambda{m}(j,:)*net.biases{m} + 0;
    intersection_low{j} =  Omega{m}(j,:) *net.biases{m} + 0;
    k = m+1;
    while k>1
        k = k-1;
        
        if k>1
            len = n{k-1};
        else
            len = n0;
        end
        
        
        
        for i = 1:len
            
            if Lambda{k}(j,:)*net.weights{k}(:,i)>0
                if k>1
                    lambda{k-1}(j,i) = param_u_linear{k-1}(1,i);    %alph_upp{k-1}(i,1);
                    Delta{k-1}(i,j)  = param_u_linear{k-1}(2,i);    %alphabet_upp{k-1}(i,1);
                    
                elseif k==1
                    lambda_0(j,i) = 1;
                end
            else
                if k>1
                    lambda{k-1}(j,i) = param_l_linear{k-1}(1,i);   %alph_low{k-1}(i,1);
                    Delta{k-1}(i,j)  = param_l_linear{k-1}(2,i);   %alphabet_low{k-1}(i,1);
                    
                elseif k==1
                    lambda_0(j,i) = 1;
                end
            end
            
            
            
            if Omega{k}(j,:)*net.weights{k}(:,i)>0
                if k>1
                    omega{k-1}(j,i) = param_l_linear{k-1}(1,i);   %alph_low{k-1}(i,1);
                    Theta{k-1}(i,j) = param_l_linear{k-1}(2,i);   %alphabet_low{k-1}(i,1);
                  
                elseif k==1
                    omega_0(j,i) = 1;
                end
            else
                if k>1
                    omega{k-1}(j,i) = param_u_linear{k-1}(1,i);   %alph_upp{k-1}(i,1);
                    Theta{k-1}(i,j) = param_u_linear{k-1}(2,i);   %alphabet_upp{k-1}(i,1);
                 
                elseif k==1
                    omega_0(j,i) = 1;
                end
            end
            
            
            
            
        end
        if k>1
            Lambda{k-1}(j,:) =  (Lambda{k}(j,:)*net.weights{k}).*lambda{k-1}(j,:);
            Omega{k-1}(j,:) =  (Omega{k}(j,:)*net.weights{k}).*omega{k-1}(j,:);
            intersection_upp{j} = intersection_upp{j}+ Lambda{k-1}(j,:)*net.biases{k-1} + Lambda{k}(j,:)*net.weights{k}*Delta{k-1}(:,j);
            intersection_low{j} = intersection_low{j}+ Omega{k-1}(j,:) *net.biases{k-1} + Omega{k}(j,:) *net.weights{k}*Theta{k-1}(:,j);
        elseif k ==1
            Lambda_0(j,:) =  (Lambda{1}(j,:)*net.weights{1}).*lambda_0(j,:);
            Omega_0(j,:) =  (Omega{1}(j,:)*net.weights{1}).*omega_0(j,:);
        end
        
    end
    
    Lam_pos = 0.5*(Lambda_0(j,:)+abs(Lambda_0(j,:)));
    Lam_neg = 0.5*(Lambda_0(j,:)-abs(Lambda_0(j,:)));
    Ubm(j,1) = Lam_pos*u0+Lam_neg*l0+intersection_upp{j};
    Ome_pos = 0.5*(Omega_0(j,:)+abs(Omega_0(j,:)));
    Ome_neg = 0.5*(Omega_0(j,:)-abs(Omega_0(j,:)));
    Lbm(j,1) = Ome_pos*l0+Ome_neg*u0+intersection_low{j};    
    
end
   
end



function [alpha_param,beta_param] =   piecewise_linear_slopes(l, u, Layer)  %%% it removes the purelin cells first and then generates the bounds


%%% This function receives the pre-activation bounds (l, u) on a single
%%% layer (Layer) and returns the bounds over the slope of activation
%%% function(alpha_param,beta_param).


if strcmp(char(Layer(1)), 'purelin')
    nln_firstloc = 0;
    dd = true;
else
    nln_firstloc = 1;
    dd = false;
end

while dd
    nln_firstloc = nln_firstloc+1;
    if ~strcmp(char(Layer(nln_firstloc)), 'purelin')
        dd = false;
    end
end
nln_locto_end = length(Layer)-nln_firstloc;
Layer = Layer(nln_firstloc:end);
l = l(end-nln_locto_end:end,1);
u = u(end-nln_locto_end:end,1);

for i = 1:length(Layer)
    if strcmp(Layer(i), 'poslin')
        alpha_param(i,1)  =  double(l(i,1)>0) ;
        beta_param(i,1)   =  double(u(i,1)>=0);
    elseif strcmp(Layer(i), 'purelin')
        clc
        disp('WARNING: you are adding purelin layers in your Lipstchitz constant computation, press any key if you want to continue ...')
        pause;
        alpha_param(i,1) = 1;
        beta_param(i,1) = 1;
    elseif strcmp(Layer(i), 'tanh')
        diffil = 1-tanh(l(i,1))^2;
        diffiu = 1-tanh(u(i,1))^2;
        diffi0 =  ((l(i,1)*u(i,1))<=0)*1;   %%% All of the activation functions have their maximum slope in zero
        alpha_param(i,1)  = min(diffil,diffiu);
        beta_param(i,1)   = max(diffi0, max(diffil,diffiu));
    elseif strcmp(Layer(i), 'sigmoid')
        diffil = sigmoid(l(i,1))* (1-sigmoid(l(i,1)));
        diffiu = sigmoid(u(i,1))* (1-sigmoid(u(i,1)));
        diffi0 =  ((l(i,1)*u(i,1))<=0)*(0.25);   %%% All of the activation functions have their maximum slope in zero
        alpha_param(i,1)  = min(diffil,diffiu);
        beta_param(i,1)   = max(diffi0, max(diffil,diffiu));
    else 
        error('the proposed activation function is not supported')
    end
end
    
end