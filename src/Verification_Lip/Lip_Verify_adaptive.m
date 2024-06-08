function [decision, Letsgo, counter_example]  =  Lip_Verify_adaptive( net, border, lb, ub, M, method_Trajectory, method_STL, LipSDP_type, Lip , status )

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


if status~=0 %% Just to meke it certain, unreliable Lip is not involved.
    Lip=[];
end

eps =  norm(ub-lb,2)/2;

decision  =  'Yes';
Letsgo = 1;
counter_example = [];

dim = size(net.weights{1},2);
for i = 1:dim
    In_axis{i} = linspace(lb(i),ub(i),M(i)+1);
end
Center = cell(M);
epsilon = cell(M);


v = ones(1,dim);
vLim = M.*v;
ready = false;

while ~ready
    
    eval(['Center{' array2index(v) '} = []']);
    eval(['epsilon{' array2index(v) '} = []']);
    for i = 1:dim
        eval(['Center{' array2index(v) '} = [Center{' array2index(v) '}; 0.5*(In_axis{' num2str(i) '}(v(' num2str(i) ')+1)+In_axis{' num2str(i) '}(v(' num2str(i) ')))];'])
        eval(['epsilon{' array2index(v) '} = [epsilon{' array2index(v) '}; 0.5*(In_axis{' num2str(i) '}(v(' num2str(i) ')+1)-In_axis{' num2str(i) '}(v(' num2str(i) ')))];'])
    end
    eval(['Center{' array2index(v) '}'])
    if status~=0
        [Lb, Ub, alpha_params, beta_params] =  Trapezius_bounds_General( eval(['Center{' array2index(v) '}']), eval(['epsilon{' array2index(v) '}']), net, border,  method_Trajectory, method_STL);
        [Lip, status, time]  =  Trapezius_Lip_SDP( net, border, alpha_params, beta_params , Lb, Ub,  LipSDP_type );
        if (status~=0)||(Lip>1000)  %%% Just to make it certain, an unreliable Lip is not involved.
            Lip=[];
            status=1;
        end
    end
    eval([ '[Rho2(' array2index(v) ') , c_e]  =  certificate_2(net, Center{' array2index(v) '}, epsilon{' array2index(v) '});'])
    
    if status==0
        certificate =  eval(['Rho2(' array2index(v) ')/Lip']);
    else
        certificate=rand;
    end
    
    if eval(['Rho2(' array2index(v) ')<0'])
        decision =  'No';
        Letsgo = 0;
        counter_example =  c_e;
    end
    
    if (status~=0)||(eps>certificate)
        [decision, Letsgo, counter_example]  =  Lip_Verify_adaptive(net, border, eval(['Center{' array2index(v) '}-epsilon{' array2index(v) '}']), eval(['Center{' array2index(v) '}+epsilon{' array2index(v) '}' ] ), M, method_Trajectory, method_STL, LipSDP_type, Lip, status );
    end
    

    if Letsgo == 0
        break;
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
end

end


function [The_min , c_e]  =  certificate_2( net, C, e )
dim = size(net.weights{1},2);
v = ones(1,dim);
vLim = 2*v;
ready = false;
j = 0;
while ~ready
    j = j+1;

    Input =  C   +   ((-1).^(v')).*e;
    rho2(j) =  NN( net, Input );


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

[The_min , index ] = min(rho2);
index=index(1)-1;
vv = dec2base(index,2);
vvv = double(vv)-48;
vvvv = vvv+1;
v = [ones(1, dim-length(vvvv)) ,  vvvv];
% disp(C)
% disp(v')
% disp(e)
% pause(2)
c_e = C + ((-1).^(v')).*e;

end




function char_indices = array2index( v )
    
    char_indices=num2str(v(1));
    for i=2:length(v)
        char_indices=[char_indices ',' num2str(v(i)) ];
    end
end