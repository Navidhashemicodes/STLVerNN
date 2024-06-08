function [Lb, Ub, alpha_params, beta_params] =  Trapezius_bounds_allReLU(Center, epsilon, border, net , method)

%%%% This function receives the set of initial states (Center, epsilon) and the relu network Trapezium-FFNN and computes the pre-activation function for all the layers, with two diffent methods.
%%% 1- exact-star : in this method the function returns the exact pre-activation bounds (Lb, Ub).
%%% 2- approx-star: in this method the function returns a relatively accurare pre-activation bound (Lb, Ub). 
%%%% The function also computes the bounds over the slope of activation functions. ( alpha_params, beta_params ) .

Net = Linear_Nonlinear(net, border);

l0 = Center-epsilon;
u0 = Center+epsilon;

STL_init = Star();
STL_init.V =  [Center eye(size(Center,1))];
STL_init.C =  zeros(1,size(Center,1));
STL_init.d =  0;
STL_init.predicate_lb = -epsilon;
STL_init.predicate_ub =  epsilon;
STL_init.dim =  size(Center,1);
STL_init.nVar =  size(Center,1);

[Lb, Ub] = star_set(STL_init, Net, method);
for i = 1:length(Net.layers)
    [alpha_params{i} ,beta_params{i}] =   piecewise_linear_slopes(Lb{i}, Ub{i}, Net.layers{i});
end

end


function [LLb, UUb] = star_set(In, Net, analysis_type)
in{1} = In;
for i = 1:length(Net.layers)

    disp([ analysis_type ' reachability analysis for Layer ' num2str(i) ' of ' num2str(length(Net.layers)-1+1) ' in Trapezius network.'])
    lenS = length(in{i});
    lenAC = 0;
    for k = 1:lenS
        in_layer{i}(k) = affineMap(in{i}(k), Net.weights{i}, Net.biases{i});
        In_layer = in_layer{i}(k);
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
        if size(In_layer.state_lb , 1 )~=0
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
                Out_layer(j).V = [in_layer{i}(k).V(1:index-1,:)  ;  Out_layer(j).V];   %%% you should check if the predicate updating process affects the V or not, I guess No, but you need to check
            elseif strcmp(analysis_type, 'approx-star')
                d_size = size(Out_layer(j).V,2)-size(in_layer{i}(k).V(1:index-1,:), 2);
                Out_layer(j).V = [[in_layer{i}(k).V(1:index-1,:) zeros(index-1, d_size) ] ;  Out_layer(j).V];
            end
            
            if size(in_layer{i}(k).state_lb , 1 )~=0
                if size(Out_layer(j).state_lb , 1 )~=0
                    Out_layer(j).state_lb = [in_layer{i}(k).state_lb(1:index-1,1) ; Out_layer(j).state_lb];
                    Out_layer(j).state_ub = [in_layer{i}(k).state_ub(1:index-1,1) ; Out_layer(j).state_ub];
                end
            else
                if size(Out_layer(j).state_lb , 1 )~=0
                    SSS =  getBox(in_layer{i}(k));
                    Out_layer(j).state_lb = [SSS.lb(1:index-1,1) ; Out_layer(j).state_lb];
                    Out_layer(j).state_ub = [SSS.ub(1:index-1,1) ; Out_layer(j).state_ub];
                end
            end

            Out_layer(j).dim = Out_layer(j).dim + (index-1);
            in{i+1}(lenAC+j) = Out_layer(j);
        end
        lenAC     = lenAC+ lenAC_now;
    end 
end

disp('pre-activation computation is in process ...')
%%%%%%  returns the pre-activation bounds over the layers  please remove
%%%%%%  the pre-activation computations for purelins later.
for i = 1:length(Net.layers)
    clear the_lb the_ub
    for j = 1:length(in_layer{i})
        if size(in_layer{i}(j).state_lb , 1 ) ==0
            SSS = getBox(in_layer{i}(j));
            the_lb(:,j) = SSS.lb;
            the_ub(:,j) = SSS.ub;
        else
            the_lb(:,j) = in_layer{i}(j).state_lb;
            the_ub(:,j) = in_layer{i}(j).state_ub;
        end
    end
    if size(the_lb,2)>1
        Lb{i} = min(the_lb')';
        Ub{i} = max(the_ub')';
    else
        Lb{i} = the_lb;
        Ub{i} = the_ub;
    end
end
LLb = Lb;
UUb = Ub;
end



function [alpha_param,beta_param] =   piecewise_linear_slopes(l, u, Layer)  %%% it removes the purelin cells first and then generates the bounds


if strcmp(Layer(1), 'purelin')
    nln_firstloc = 1;
    dd = true;
    while dd
        nln_firstloc = nln_firstloc+1;
        if ~strcmp(Layer(nln_firstloc), 'purelin')
            dd = false;
        end
    end
else
    nln_firstloc = 1;
end
% if strcmp(char(Layer(1)), 'purelin')
%     nln_firstloc=0;
%     dd=true;
% else
%     nln_firstloc=1;
%     dd=false;
% end
% 
% while dd
%     nln_firstloc=nln_firstloc+1;
%     if ~strcmp(char(Layer(nln_firstloc)), 'purelin')
%         dd=false;
%     end
% end

Layer = Layer(nln_firstloc:end);
l = l(nln_firstloc:end,1);
u = u(nln_firstloc:end,1);

for i = 1:length(Layer)
    if strcmp(Layer(i), 'poslin')
%         disp('ding')
        alpha_param(i,1)  =  double(l(i,1)>0) ;
        beta_param(i,1)   =  double(u(i,1)>=0);
    elseif strcmp(Layer(i), 'purelin')
        clc
        error('WARNING: you are adding purelin layers in your Lipstchitz constant computation, press any key if you want to continue ...')
        pause;
        alpha_param(i,1) = 1;
        beta_param(i,1) = 1;
    end 
end
end