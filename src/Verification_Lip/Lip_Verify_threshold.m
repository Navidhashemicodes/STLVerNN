function Loc_Lip  =  Lip_Verify_threshold( net, model_dim, cntrl_dim, border, lb, ub, method_Trajectory, method_STL, LipSDP_type , threshold )

% decision: is the main result. if decision='Yes' then the controller is verified and once decision='No' the the controller is rejected and a counter-example is returned.
% Letsgo: Once we found the counter example this variable is set to 0, which implies the process is to immidiately terminated.
% counter_example: This is the counter example we find, it is [] if the controller is verified and is nonempty if the controller is rejected.
% net: Represents the Trapezium-FFNN.
% border: Is the position for the first layer of STL FFNN on Trapezium FFNN.  
% lb: is the lower bound on the set of initial states.
% ub: is the upper bound on the set of initial states.
% M : Is an array containing the number of partitions for every dimension of the set of initial states.  
% method_Trajectory: Is the method we adopt to compute the pre-activation bound on the trajectory, The options for this method provided in function
% Trapezius_bounds_General()  are as follows:
%  1-method_Trajectory='Linear-CROWN';
%  2-method_Trajectory='Quadratic-CROWN';
%  3-method_Trajectory='Quadratic-Gurobi';
% method_STL : Is the method we apply to compute pre-activation bounds on the STL-FFNN. This method are as follows:
%  1- method_STL='approx-star' ;
%  2- method_STL='exact-star' ;
% LipSDP_type: Is the method we adopt to charactrize the convex optimization over the layers of Trapezium FFNN this methods are sorted as follows:
%  1- LipSDP_type= 'LipSDP-layer';  Is a trade off between Running time and accuracy. 
%  2- LipSDP_type= 'LipSDP-neuron'; Is the most accurate convex programming technique but sometimes suffers from Running time.  
%  3- LipSDP_type= 'LipSDP-network';Is the most efficient convex programming technique but sometimes sufferes from accuracy, and even worse, it may lead to infeasibility. 


dim = size(net.weights{1},2);
% depth=length(net.layers);
% portion_size=threshold/depth;
% Num_disc=zeros(1,dim);
In_axis=cell(1,dim);
% Total=ub-lb;
for i=1:dim
%     Num_disc(i)=ceil(Total(i)/portion_size);
%     In_axis{i}=linspace(lb(i),ub(i), Num_disc(i)+1);
    In_axis{i}=linspace(lb(i),ub(i), threshold(i)+1);
end

sorter=0;
v = ones(1,dim);
vLim = threshold;
ready = false;
while ~ready
    sorter=sorter+1;
    eval(['Center{' array2index(v) '} = []']);
    eval(['epsilon{' array2index(v) '} = []']);
    for i = 1:dim
        if length(In_axis{i})>=2
            eval(['Center{' array2index(v) '} = [Center{' array2index(v) '}; 0.5*(In_axis{' num2str(i) '}(v(' num2str(i) ')+1)+In_axis{' num2str(i) '}(v(' num2str(i) ')))];'])
            eval(['epsilon{' array2index(v) '} = [epsilon{' array2index(v) '}; 0.5*(In_axis{' num2str(i) '}(v(' num2str(i) ')+1)-In_axis{' num2str(i) '}(v(' num2str(i) ')))];'])
        else
            eval(['Center{' array2index(v) '} = [Center{' array2index(v) '}; In_axis{' num2str(i) '}(v(' num2str(i) '))];'])
            eval(['epsilon{' array2index(v) '} = [epsilon{' array2index(v) '}; 0];'])
        end
    end
    [Lb, Ub, alpha_params, beta_params] =  Trapezius_bounds_General( eval(['Center{' array2index(v) '}']), eval(['epsilon{' array2index(v) '}']), net, model_dim, cntrl_dim, border,  method_Trajectory, method_STL);
    eval(['[Lip(' array2index(v) '), status(' array2index(v) '), time(' array2index(v) ')]  =  Trapezius_Lip_SDP( net, border, alpha_params, beta_params , Lb, Ub,  LipSDP_type );'])
    eval(['LIP(sorter)=Lip(' array2index(v) ');' ]);
    eval(['STATUS=status(' array2index(v) ');' ]);
    if STATUS ~=0
        error([num2str(threshold), ' was not enough, try a smaller threshold (updated default is currently 12.5).'])
    end
    
    
    ready = true;
    for k = 1:dim                                   %index updater
        v(k) = v(k) + 1;
        if v(k) <=  vLim(k)
            ready = false;
            break;
        end
        v(k) = 1;
    end
    
end

Loc_Lip=max(LIP);

end



function char_indices = array2index( v )
    
    char_indices=num2str(v(1));
    for i=2:length(v)
        char_indices=[char_indices ',' num2str(v(i)) ];
    end
end