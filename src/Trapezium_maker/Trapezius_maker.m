function  [Net, border]  = Trapezius_maker(model, nncontroller, predicate_transformation , str)
%%%This function receives the model (FFNN or LTV) and controller (FFNN) and produces a single FFNN that represents all the trajectory. 
%%%The length of the trajectory is also determined from the STL specifications. In the next attempt The function generates a FFNN from STL
%%%specifications (STL-FFNN) and connects STL-FFNN to the trajectory by means of bariers for predicates introduced to the function. The result of the 
%%% function is the network --->Net<---- Which receives the initial states,
%%% simulates the trajectory and resturns its robustness.  The out put
%%% --->border<---- is the index of the layer where STL-FFNN starts.

[s1,s2] = size(str);
TT = regexprep(strjoin(reshape(string(str)', [1,s1*s2]))  , ' ',  ''); %%% If you print TT you think it is the STL formula written in terms of '&' and '|'.
Tree = STL2Tree(char(TT));  %%%Generates the Tree from the STL specifications which is introduced as char array.

net = Tree2FFNN(Tree);  %%% Generates the STL FFNN from the tree. 


%%%%% modification of the net for repeatative inputs
%%%%% The logic Tree can contains repeatative leaf-nodes The following lines 
%%%%% Detects the repeatative leaf nodes and removes them, it also modifys
%%%%% the first weight matrix of STL-FFNN with merging the rows
%%%%% corresponding to removed lead-nodes.

n = nnodes(Tree);  %%% This function returns the number of nodes on the Tree.
ids = depthfirstiterator(Tree, 1, false);  %%% This function sorts the nodes of the Tree based on DFS traversal.
leaf_IDs = ids;  %%% Here we want to draw the leaf nodes. from the nodes. The leaf nodes are successfully drawn out after the while loop. 
i = 0;
while i<n
    i = i+1;
    if ~isleaf(Tree,leaf_IDs(i))
        leaf_IDs(i) = [];
        i = i-1;
        n = n-1;
    end
end
leaf_values = cell(1,length(leaf_IDs));   %%% when the lead nodes are extracted we return the value of the nodes and sort 
                                          %%%them in the cell leaf_values. The following for loop is taking care of this step.
for i = 1:length(leaf_IDs)
    leaf_values{i} = {get(Tree, leaf_IDs(i))};
end

W1 = sparse(net.weights{1});
i = 0;
while i<length(leaf_values)   %%%% We have repeatative leaf nodes on the tree, what I mean is, there exist some predicates who are contributing in diffenet formulas. We 
                              %%%% need to collect them in one input in FFNN and consequently their wights should also be merged, the following while loop takes care of this process.
    i = i+1;
    j = i;
    while j<length(leaf_values)
        j = j+1;
        if strcmp(leaf_values{i},leaf_values{j})
            leaf_values(j) = [];
%             leaf_values=[leaf_values(1:j-1), leaf_values(j+1:end)];
            W1(:,i) = sparse(W1(:,i)+W1(:,j));
            W1(:,j) = [];
            j = j-1;
        end
    end
end


%%%% Extract the time horizon: In the followinf for loop we extract the
%%%% time horizon from the STL specification.
ss=0;
for i = 1:length(leaf_values)
    character = char(leaf_values{i});
    tt = str2num(character);
    ss = max(ss, tt(2));
end
T = ss;


net.weights{1} = sparse(W1);   %%% as you see net represents STL-FFNN and Net represent the Trapezium-FFNN. Here the first weight of Net is successfully updated for merging repeatative nodes.  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch model(1).type   %%%% Now that STL-FFNN is finished we focus on the Trajectory part of Trapezium-FFNN.
    
    
    
    
    
    
    
    case 'LTV'
        %%%%  NN controller accepts O{t}=C{t}s{t}+D{t} as input  a=FFNN(O),
        %%%%  s(t+1)=A{t}s(t) + B{t} a(t) + E{t}
        C = model.C;
        D = model.D;  
        A = model.A;
        B = model.B;
        E = model.E;
        
        system_dim = size(A{1},1);
        
        ell = length(nncontroller(1).layers);
        weights = nncontroller(1).weights;
        zp_end_size = size(weights{ell+1},2);
        
        
        
        for k = 1:T  %%% Given the FFNN for controller and the Linear model we can construct a FFNN that represents all the trajectory. Consider the states should 
                   %%% also be transferred once computed upto the last layer. This implies we need to add also linear activations. the following loop performs this operation. 
                   
                   
            weights = nncontroller(k).weights;
            biases  = nncontroller(k).biases;
            layers  = nncontroller(k).layers;
            
            
            index=(k-1)*ell;
            if k==1
                Net.weights{index+1} = sparse([eye(system_dim) ; weights{1}*C{k}]);
                Net.biases{index+1} = sparse([zeros(system_dim,1) ; weights{1}*D{k} + biases{1}]);
                
            else
                Net.weights{index+1} = sparse(blkdiag(eye((k-2)*system_dim) , [eye(system_dim)       zeros(system_dim ,  zp_end_size)      ;...
                                                                                    A{k}                    B{k}*weights{ell+1}           ;...
                                                                            weights{1}*C{k}*A{k}   weights{1}*C{k}*B{k}*weights{ell+1}    ] ));

                Net.biases{index+1} = sparse([zeros((k-1)*system_dim , 1); B{k}*biases{ell+1}+E{k} ;...
                                              weights{1}*C{k}*(B{k}*biases{ell+1}+E{k})+weights{1}*D{k}+biases{1}]);
            end
            clear LL
            LL = cell(k*system_dim+size(biases{1},1),1);  %%% our structure to present a FFNN is different from the public structures, we present weights, biases and we also 
                                                         %%% represent the layers with a char array, that contains the name of activations, we further emply str2func() to make
                                                         %%% the function from this characters. the parameter LL is introduced for this purpose
            LL(1:k*system_dim) = {'purelin'};
            LL(k*system_dim+1:end) = layers{1};
            Net.layers{index+1} = LL ;
            for j = 2:ell
                Net.weights{index+j} = sparse(blkdiag(eye(k*system_dim) , weights{j}));
                Net.biases{index+j} = sparse([zeros(k*system_dim,1) ; biases{j} ]) ;
                
                clear LL
                LL = cell(k*system_dim+size(biases{j},1),1);
                LL(1:k*system_dim) = {'purelin'};
                LL(k*system_dim+1:end) = layers{j};
                Net.layers{index+j} = LL ;
            end
        end
        index = T*ell;
        Net.weights{index+1} = sparse(blkdiag(eye((T-1)*system_dim) , [eye(system_dim)  zeros(system_dim , zp_end_size)     ;...
                                                                           A{k}              B{k}*weights{ell+1}          ]  ));
                                                                 
        Net.biases{index+1} = sparse([  zeros(T*system_dim,1)    ;   B{k}*biases{ell+1}+E{k}  ]) ;        
        
        
        
        index = T*ell;
        border = index+1;   %%% border is the layer where the STL-FFNN starts on Trapezium-FFNN.
        
        
        
        
        
        
    case 'FFNN'   %%% Given the FFNN for controller and the FFNN for model we can construct a FFNN that represents all the trajectory. Consider the states should 
                  %%% also be transferred once computed upto the last layer. This implies we need to add also linear activations. the following loop performs this operation.
        
        %%%%  NN controller accepts O{t}=C{t}s{t}+D{t} as input  a=FFNN(O),
        
        system_dim = size(model(1).C,2); 
        ell = length(nncontroller(1).layers);
        
        
        
        ell_m = length(model(1).layers);
        weights = nncontroller(1).weights;
        zp_end_size_m = size(weights{ell+1},2);
        
        
        for k = 1:T
            
            weights = nncontroller(k).weights;
            biases  = nncontroller(k).biases;
            layers  = nncontroller(k).layers;
            
            
            weights_m = model(k).weights;
            biases_m  = model(k).biases;
            layers_m  = model(k).layers;
            
            
            index = (k-1)*(ell_m+ell);
            if k==1
                Net.weights{index+1} = sparse([eye(system_dim) ; weights{1}*model(k).C]);
                Net.biases{index+1}  = sparse([zeros(system_dim,1) ; weights{1}*model(k).D + biases{1}]);
            else
                Net.weights{index+1} = sparse(blkdiag(eye((k-1)*system_dim) , [weights_m{end} ; weights{1}*model(k).C*weights_m{end}]));
                Net.biases{index+1}  = sparse([zeros((k-1)*system_dim ,1)  ;  biases_m{end}   ;  weights{1}*model(k).C*biases_m{end}+weights{1}*model(k).D+biases{1} ]);
            end
            clear LL
            LL = cell(k*system_dim+size(biases{1},1),1);   %%% our structure to present a FFNN is different from the public structures, we present weights, biases and we also 
                                                          %%% represent the layers with a char array, that contains the name of activations, we further emply str2func() to make
                                                          %%% the function from this characters. the parameter LL is introduced for this purpose
            LL(1:k*system_dim) = {'purelin'};
            LL(k*system_dim+1:end) = layers{1};
            Net.layers{index+1} = LL ;
            
            for j = 2:ell
                Net.weights{index+j} = sparse(blkdiag(eye(k*system_dim) , weights{j}));
                Net.biases{index+j} = sparse([zeros(k*system_dim,1) ; biases{j} ]) ;
                
                clear LL
                LL = cell(k*system_dim+size(biases{j},1),1);
                LL(1:k*system_dim) = {'purelin'};
                LL(k*system_dim+1:end) = layers{j};
                Net.layers{index+j} = LL ;
                
            end
            
            index = index+ell;
            Net.weights{index+1} = sparse(blkdiag(eye((k-1)*system_dim)   ,   [    eye(system_dim)          zeros(system_dim,zp_end_size_m)   ;... 
                                                                                    weights_m{1}*blkdiag(eye(system_dim),weights{end})      ]  ));
                                                                                                                                   
            Net.biases{index+1} = sparse([zeros(k*system_dim,1)  ;  weights_m{1}*[zeros(system_dim,1);biases{end}]+biases_m{1}]);
            
            clear LL
            LL = cell(k*system_dim+size(biases_m{1},1),1);
            LL(1:k*system_dim) = {'purelin'};
            LL(k*system_dim+1:end) = layers_m{1};
            Net.layers{index+1} = LL ;
            
            
            for j = 2:ell_m
                Net.weights{index+j} = sparse(blkdiag(eye(k*system_dim) , weights_m{j}));
                Net.biases{index+j} = sparse([zeros(k*system_dim,1) ; biases_m{j} ]) ;
                
                clear LL
                LL = cell(k*system_dim+size(biases_m{j},1),1);
                LL(1:k*system_dim) = {'purelin'};
                LL(k*system_dim+1:end) = layers_m{j};
                Net.layers{index+j} = LL ;
            end
        end
        index = T*(ell+ell_m);
        Net.weights{index+1} = sparse(blkdiag(eye(T*system_dim) , weights_m{end}));
        Net.biases{index+1} = sparse([zeros(k*system_dim,1) ; biases_m{end} ]) ;
        
        
        
        index = T*(ell+ell_m);
        border = index+1;
        
          
        
        
        
        
        
        
        

    otherwise 
        error('we only support Linear models or Feed forward neural network models')
end






Transformation = zeros(length(leaf_values), (T+1)*system_dim);
threshold = zeros(length(leaf_values), 1);
for i = 1:length(leaf_values)
    character = char(leaf_values{i});    
    tt = str2num(character);            
    what_predicate = tt(1);         %%%%  Takes the element of leaf_value, this represent the number of predicate    
    what_time_State = tt(2);        %%%%  Takes the time that predicate is defined on.
    time_index = what_time_State*system_dim;
    if what_time_State>T+1 %%% This if is not necessary. It wont happen at all , we can remove it.
        error('there is no state computed for this predicate.')
    end
    Transformation(i,time_index+1:time_index+system_dim) = predicate_transformation.C(what_predicate,:);
    threshold(i,1)                                       = predicate_transformation.d(what_predicate,1);
end

%%% map the states to leaf predicates and then on the first layer of STL NN

Net.weights{index+1} = sparse(net.weights{1}*Transformation*Net.weights{index+1});  %%% This is the linear map for predicates that that connects the Trsjectory to the STL-FFNN
                                                                                  %%% as you see we didnt simply multiply last weight of trajectory in the first weight of the STL-FFNN. 
Net.biases{index+1} = sparse(net.weights{1}*Transformation*Net.biases{index+1} + net.weights{1}*threshold);   %%% The linear map also tells us the bias vector needs to be update like this
Net.layers{index+1} = net.layers{1};    %%%% The activation functions does not change. Only the weights and biases will change in connection process.

for i=2:length(net.layers)
    Net.layers{index+i} = net.layers{i};
    Net.weights{index+i} = sparse(net.weights{i});
    Net.biases{index+i} = sparse(zeros(size(net.weights{i},1) , 1 ));
end
index = index+length(net.layers);
Net.weights{index+1} = sparse(net.weights{end});
Net.biases{index+1} = sparse(zeros(size(net.weights{end},1) , 1 ));  %%% As there exist no bias vector on the STL-FFNN so we assign 0 to their bias verctors.

end