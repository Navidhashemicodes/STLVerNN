function  interval = Trapezius_ReLUplex_verifier_withNet(Center, epsilon, Net, analysis_type, num_Cores)

%%%% This code is in fact   Trapezius_ReLUplex_verifier.m  but I made it
%%%% such that it doesnt compute the Trapezium-FFNN anymore. This makes the
%%%% code faster to be run. There is no need to check this code as it is in
%%%% fact Trapezius_ReLUplex_verifier.m



system_dim = size(Center,1);

init_star_state = Star();     %%% consider the following lines, this how I initiate the Star() set, my experiments shows, if I initiate Star() set in another way, then I will face computational errors. 
init_star_state.V =  [Center eye(system_dim)];
init_star_state.C =  zeros(1,system_dim);
init_star_state.d =  0;
init_star_state.predicate_lb = -epsilon;
init_star_state.predicate_ub =  epsilon;
init_star_state.dim =  system_dim;
init_star_state.nVar =  system_dim;

in{1} = init_star_state;   %%%% This is the Star() set that will be introduced to the neural network. I am thinking to upgrade the in away it accept multiple Star() sets. In that case we can also include non-convex sets in our research.
check_bound = true;    %%%%
len = length(Net.layers);
for i = 1:len
    disp(['analysis of Layer ' num2str(i) ' of ' num2str(len) '.'])
    lenS = length(in{i});
    lenAC = 0;
    
    for k = 1:lenS
        if lenS>250
            interval = 'Unknown';
            check_bound = false;
            break
        else
        
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
        %%%%    'exact-star'   or    'approx-star'
        [Out_layer, ~] = reach(F,In_layer, analysis_type, num_Cores);
        
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
        
            check_bound = true;
        end
    end
    if check_bound==false
        break;
    end
end


%%%%%%  returns the pre-activation bounds over the layers  please remove
%%%%%%  the pre-activation computations for purelins later.
% for i=1:len
%     clear the_lb the_ub
%     for j=1:length(in_layer{i})
%         if size(in_layer{i}(j).state_lb , 1 ) ==0
%             SSS = getBox(in_layer{i}(j));
%             the_lb(:,j)=SSS.lb;
%             the_ub(:,j)=SSS.ub;
%         else
%             the_lb(:,j) = in_layer{i}(j).state_lb;
%             the_ub(:,j) = in_layer{i}(j).state_ub;
%         end
%     end
%     if size(the_lb,2)>1
%         Lb(:,i)=min(the_lb')';
%         Ub(:,i)=max(the_ub')';
%     else
%         Lb(:,i)=the_lb;
%         Ub(:,i)=the_ub;
%     end
% end

if check_bound
lenS = length(in{len+1});
for k = 1:lenS
    range(k) = affineMap(in{len+1}(k), Net.weights{end}, Net.biases{end});
    
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
end


% figure
% for i=1:system_
%     
%     subplot(1,6,i)
%     plot(Lb(i,:));
%     hold on
%     plot(Ub(i,:));
%     
% end

