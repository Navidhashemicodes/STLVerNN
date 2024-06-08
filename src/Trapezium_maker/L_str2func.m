function post_activation = L_str2func( char_array , pre_activation)


%%%% As we mentioned in function NN.m  Our data structure to represent NN
%%%% is different from the others. we propose a char array which is the
%%%% layer. To compute the post-activation of the layer from the
%%%% pre-activations we need to change the layer's strings to the functions. This function performs this task and computes the post activation given a pre-activation. 

n = size(pre_activation,1);

post_activation = zeros(n,1);
for i = 1:n
    f = str2func(char_array{i});
    post_activation(i,1) =  f(pre_activation(i));
end


end
    