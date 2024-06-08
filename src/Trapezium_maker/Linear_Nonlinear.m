function Net=Linear_Nonlinear(net, border)

%%%% This function is written for Trapezium-FFNN. We know the location of linear and
%%%% non-linear activations are separated. Thus we dont process the layers over the trajectory 
%%%% but when we step forward to STL-FFNN they are mixed. We need to separate the location of them
%%%% This function takes the input 'border' to see where the STL-FFNN starts. then plays with weights
%%%% of STL-FFNN to separate the lication of linear and nonlinear activations. Thus the result is a NN that the linear activations are sorted on top in every single layer. 


Net = net;
len = length(net.layers);

W1 = net.weights{border};
for i = border:len
    W2 = net.weights{i+1};
    b = net.biases{i};
    L = net.layers{i};
    ss = 0;
    pure_location = [];
    W1_new = W1;
    W2_new = W2;
    b_new = b;
    L_new = L;
    for j = 1:length(L)
        if strcmp(L(j), 'purelin')
           ss = ss+1; 
           pure_location = [pure_location j];
           W1_new(ss,:) = W1(j,:);
           b_new(ss,1) = b(j,1);
           L_new(ss) = L(j);
           W2_new(:,ss) = W2(:,j);
        end
    end
    TT = 1:length(L);
    nonlin_location = setdiff(TT,pure_location);
    index = length(pure_location);
    for j = 1:length(nonlin_location)
        W1_new(index+j,:) = W1(nonlin_location(j), :);
        b_new(index+j,1)  = b(nonlin_location(j),1)  ;
        L_new(index+j)    = L(nonlin_location(j))    ;
        W2_new(:,index+j) = W2(:, nonlin_location(j));
    end
    Net.weights{i} = W1_new;
    Net.biases{i} = b_new;
    Net.layers{i} = L_new;
    W1 = W2_new;
end
Net.weights{len+1} = W2_new;
Net.biases{len+1} = net.biases{len+1};



end