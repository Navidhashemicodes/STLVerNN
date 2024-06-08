function y=NN(Net,x)
%%%% As you see in this toolbox we present our own data structure for NN. We return a NN as a tuple (weights, biases, layers).
%%%%  1-weights: is a cell of weight elements over the NN.
%%%%  2- bises : is a cell of bias vector elements in  NN.
%%%%  3- layers: is char array, which contains the name of activation functions. we can easily convert this char arrays to function using str2fuc().
out = x;
for i = 1:length(Net.layers)
    out = L_str2func(Net.layers{i},Net.weights{i}*out+Net.biases{i});
end
y = Net.weights{end}*out+Net.biases{end};


end