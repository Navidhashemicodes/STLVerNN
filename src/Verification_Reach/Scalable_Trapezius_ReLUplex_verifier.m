function  interval = Scalable_Trapezius_ReLUplex_verifier(Center, epsilon, model, nncontroller, predicate_transformation, str, analysis_type, pred_thresh, num_thresh, neuro_thresh)

%%% This function is in fact Trapezius_ReLUplex_verifier.m but there is three new inputs which will be discussed here.
%%%   1- pred_thresh  :  We monitor the evolution of star sets. If the number of the predicates exceeds prde_thresh we replace that star set with a box .
%%%   2- num_thresh   :  We monitor the evolution of star sets. If the number of the star sets exceeds num_thresh we replace the union of them with a box .
%%%   3- neuro_thresh :  If the number of relu activation functions on a specific layer, exceeds neur_thresh, then we split the layer in layes with size atmost neuro_thresh. We then pass the star sets through the splitted layers.
 
system_dim = size(Center,1);

[net, border] = Trapezius_maker(model, nncontroller, predicate_transformation , str);
Net = Linear_Nonlinear(net, border);

init_star_state = Star();
init_star_state.V =  [Center eye(system_dim)];
init_star_state.C =  zeros(1,system_dim);
init_star_state.d =  0;
init_star_state.predicate_lb = -epsilon;
init_star_state.predicate_ub =  epsilon;
init_star_state.dim =  system_dim;
init_star_state.nVar =  system_dim;

in{1} = init_star_state;

len = length(Net.layers);

for i = 1:len
    disp(['analysis of Layer ' num2str(i) ' of ' num2str(len) '.'])
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
    %%%%index is the position of the first nonlinear
    the_length = length(Net.layers{i})-index+1;
    
    if the_length>neuro_thresh
        jj = 1;
        lenS = length(in{i});
        for k = 1:lenS
            Inn_l(k) = affineMap(in{i}(k), full(Net.weights{i}), full(Net.biases{i}));
        end
        lenS_u = lenS;
        while jj*neuro_thresh <=  the_length
            for k = 1:lenS_u
                Inn_l_u(k).V = Inn_l(k).V(index:index+neuro_thresh-1, :);   %%% index is always the first position for the nonlinear for all segments
                if size(Inn_l(k).state_lb , 1 )~=0
                    Inn_l_u(k).state_lb = Inn_l(k).state_lb(index:index+neuro_thresh-1, 1);
                    Inn_l_u(k).state_ub = Inn_l(k).state_ub(index:index+neuro_thresh-1, 1);
                end
                Inn_l_u(k).dim =  size(Inn_l_u(k).V , 1 );
            end
            
            lenAC = 0;
            for k = 1:lenS_u
                
                Layers = LayerS(eye(Inn_l_u(k).dim), zeros(Inn_l_u(k).dim,1) , 'poslin');
                Ln     = LayerS(eye(Inn_l_u(k).dim), zeros(Inn_l_u(k).dim,1) , 'purelin');
                Layers = [Layers Ln];
                
                F = FFNNS(Layers);
                
                %%%%    'exact-star'   or    'approx-star'
                [Out_layer, ~] = F.reach(Inn_l_u(k), analysis_type);
                
                lenAC_now = length(Out_layer);
                for j = 1:lenAC_now
                    if strcmp(analysis_type, 'exact-star')
                        Out_layer(j).V = [Inn_l(k).V(1:index-1,:)  ;  Out_layer(j).V ; Inn_l(k).V(index+neuro_thresh:end,:)];   %%% you should check if the predicate updating process affects the V or not, I guess No, but you need to check
                    elseif strcmp(analysis_type, 'approx-star')
                        d_size = size(Out_layer(j).V,2)-size(Inn_l(k).V, 2);
                        size2 =  size(Inn_l(k).V(index+neuro_thresh:end,:) ,1);
                        Out_layer(j).V = [[Inn_l(k).V(1:index-1,:) zeros(index-1, d_size) ] ;  Out_layer(j).V  ;   [ Inn_l(k).V(index+neuro_thresh:end,:)    zeros(size2, d_size)  ]  ];
                    end
                    
                    if size(Inn_l(k).state_lb , 1 )~=0
                        if size(Out_layer(j).state_lb , 1 )~=0
                            Out_layer(j).state_lb = [Inn_l(k).state_lb(1:index-1,1) ; Out_layer(j).state_lb ; Inn_l(k).state_lb(index+neuro_thresh:end,1)];
                            Out_layer(j).state_ub = [Inn_l(k).state_ub(1:index-1,1) ; Out_layer(j).state_ub ; Inn_l(k).state_ub(index+neuro_thresh:end,1)];
                        end
                    else
                        if size(Out_layer(j).state_lb , 1 )~=0
                            SSS =  getBox(Inn_l(k));
                            Out_layer(j).state_lb = [SSS.lb(1:index-1,1) ; Out_layer(j).state_lb ; SSS.lb(index+neuro_thresh:end,1)];
                            Out_layer(j).state_ub = [SSS.ub(1:index-1,1) ; Out_layer(j).state_ub ; SSS.ub(index+neuro_thresh:end,1)];
                        end
                    end
                    
                    Out_layer(j).dim = Inn_l(k).dim;
                    if Out_layer(j).nVar > pred_thresh
                        Out_layer(j) = rectangle_bounder(Out_layer(j));
                    end
                    Inn_l_next(lenAC+j) = Out_layer(j);
                end
                lenAC     = lenAC+ lenAC_now;
            end
            
            
            lenS_u = length(Inn_l_next);
            
            if lenS_u>num_thresh
                Inn_l_next = rectangle_bounder(Inn_l_next);
                lenS_u = 1;
            end
            Inn_l = Inn_l_next;
            jj = jj+1;
            index = index+neuro_thresh;
        end
        
        if jj*neuro_thresh~= the_length
            
            index = index-neuro_thresh;
            for k = 1:lenS_u
                Inn_l_u(k).V = Inn_l(k).V(index:end, :);
                if size(Inn_l(k).state_lb , 1 )~= 0
                    Inn_l_u(k).state_lb = Inn_l(k).state_lb(index:end, 1);
                    Inn_l_u(k).state_ub = Inn_l(k).state_ub(index:end, 1);
                end
                Inn_l_u(k).dim =  size(Inn_l_u(k).V , 1 );
            end
            
            lenAC = 0;
            
            for k = 1:lenS_u
                
                Layers = LayerS(eye(Inn_l_u(k).dim), zeros(Inn_l_u(k).dim,1) , 'poslin');
                Ln     = LayerS(eye(Inn_l_u(k).dim), zeros(Inn_l_u(k).dim,1) , 'purelin');
                Layers = [Layers Ln];
                
                F = FFNNS(Layers);
                
                %%%%    'exact-star'   or    'approx-star'
                [Out_layer, ~] = F.reach(Inn_l_u(k), analysis_type);
                
                lenAC_now = length(Out_layer);
                for j = 1:lenAC_now
                    if strcmp(analysis_type, 'exact-star')
                        Out_layer(j).V = [Inn_l(k).V(1:index-1,:)  ;  Out_layer(j).V ];   %%% you should check if the predicate updating process affects the V or not, I guess No, but you need to check
                    elseif strcmp(analysis_type, 'approx-star')
                        d_size = size(Out_layer(j).V,2)-size(Inn_l(k).V, 2);
                        Out_layer(j).V = [[Inn_l(k).V(1:index-1,:) zeros(index-1, d_size) ] ;  Out_layer(j).V  ];
                    end
                    
                    if size(Inn_l(k).state_lb , 1 )~=0
                        if size(Out_layer(j).state_lb , 1 )~=0
                            Out_layer(j).state_lb = [Inn_l(k).state_lb(1:index-1,1) ; Out_layer(j).state_lb ];
                            Out_layer(j).state_ub = [Inn_l(k).state_ub(1:index-1,1) ; Out_layer(j).state_ub ];
                        end
                    else
                        if size(Out_layer(j).state_lb , 1 )~=0
                            SSS =  getBox(Inn_l(k));
                            Out_layer(j).state_lb = [SSS.lb(1:index-1,1) ; Out_layer(j).state_lb ];
                            Out_layer(j).state_ub = [SSS.ub(1:index-1,1) ; Out_layer(j).state_ub ];
                        end
                    end
                    
                    Out_layer(j).dim = Inn_l(k).dim;
                    if Out_layer(j).nVar > pred_thresh
                        Out_layer(j) = rectangle_bounder(Out_layer(j));
                    end
                    Inn_l_next(lenAC+j) = Out_layer(j);
                end
                lenAC     = lenAC+ lenAC_now;
            end
            
            lenS_u = length(Inn_l_next);
            
            if lenS_u>num_thresh
                Inn_l_next = rectangle_bounder(Inn_l_next);
            end
            Inn_l = Inn_l_next;
        end
        in{i+1} = Inn_l;
        
    else
        lenS = length(in{i});
        lenAC = 0;
        for k = 1:lenS
            in_layer{i}(k) = affineMap(in{i}(k), full(Net.weights{i}), full(Net.biases{i}));
            
            In_layer = in_layer{i}(k);
            
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
            
            %%%%    'exact-star'   or    'approx-star'
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
                if Out_layer(j).nVar > pred_thresh
                    Out_layer(j) = rectangle_bounder(Out_layer(j));
                end
                in{i+1}(lenAC+j) = Out_layer(j);
            end
            lenAC     = lenAC+ lenAC_now;
        end
        lenS = length(in{i+1});
        
        if lenS>num_thresh
            in{i+1} = rectangle_bounder(in{i+1});
        end
        
    end
    
    
    
end


lenS = length(in{len+1});
for k = 1:lenS
    range(k) = affineMap(in{len+1}(k), full(Net.weights{end}), full(Net.biases{end}));
    
    if size(range(k).state_lb , 1 ) ==0
        SSS = getBox(range(k));
        the_lb(:,k) = SSS.lb;
        the_ub(:,k) = SSS.ub;
    else
        the_lb(:,k) = range(k).state_lb;
        the_ub(:,k) = range(k).state_ub;
    end
end

if size(the_lb,2)>1
    Lb = min(the_lb')';
    Ub = max(the_ub')';
else
    Lb = the_lb;
    Ub = the_ub;
end

interval = [ Lb , Ub ];


end



function S = rectangle_bounder(stars)
            % @stars: an array of stars
            % @S: a new star which is a convex hull of the input stars
            
            % author: Dung Tran
            % date: 1/21/2019
            
            n = length(stars);
            dim = stars(1).dim; 
            
            
            for i = 2:n
                if ~isa(stars(i), 'Star')
                    error('The %d th object is not a star', i);
                end
                if stars(i).dim ~= dim
                    error('Inconsistent dimensions between stars');
                end
            end
            
            for i = 1:n
                X(i) = stars(i).toPolyhedron;
            end
            
            U = PolyUnion(X);
            S = U.convexHull;
% S=stars.get_convex_hull(star);
% S=get_Convex_hull(I)

end
% 
% function I = rectangle_bounder(J)
% disp('ding')
% % I=get_convex_hull(J);
% I=getZono(J);
% % len=length(J);
% % 
% % for j=1:len
% %     if size(J(j).state_lb , 1 ) ==0
% %         SSS = getBox(J(j));
% %         the_lb(:,j)=SSS.lb;
% %         the_ub(:,j)=SSS.ub;
% %     else
% %         the_lb(:,j) = J(j).state_lb;
% %         the_ub(:,j) = J(j).state_ub;
% %     end
% % end
% % if size(the_lb,2)>1
% %     Lb=min(the_lb')';
% %     Ub=max(the_ub')';
% % else
% %     Lb=the_lb;
% %     Ub=the_ub;
% % end
% % 
% % system_dim=size(Lb,1);
% % Center = 0.5*(Lb+Ub);
% % epsilon= 0.5*(Ub-Lb);
% % 
% % I=Star();
% % I.V= [Center eye(system_dim)];
% % I.C= zeros(1,system_dim);
% % I.d= 0;
% % I.predicate_lb=-epsilon;
% % I.predicate_ub= epsilon;
% % I.dim= system_dim;
% % I.nVar= system_dim;
% 
% end