function [Weights, layers] = pure_poslin(type, num_inputs)

%%%%   This function gets the number of inputs who are going to contribute in an max or min process and returns the NN (weights,layers) who returns the min or the max. Consider the parameter 'type' can be 'max' or 'min' 
%%%%   First of all, I wanna make it clear, why this function is recursive,
%%%%   Assume you have written a function that do this for 2 inputs. In addition we have another function that returns a neural network for 3 inputs. what about 4? for 4 inputs we pair the inputs into piars of 2 and take 
%%%%   the 'min/max' process for each pair. So clearly in the next step we need to process 'min/max' for 2 elements which we had a function before. What about 5? Again we pair all of these 5 inputs and take 'min/max' of each pair.
%%%%   Then the problem will be reduced to take a 'min/max' from 3 elements, which we had a function for this before. We can continue this process for an arbitrary number of inputs and the problem is reducible to find the 'min/max'
%%%%   of 3 or 2. The following lines of codes implements this idea and generates NN for 'min/max' process in a recursive way.



mapper = [1 1; -1 -1 ; 1 -1 ; -1 1];
min_mapper = [0.5 -0.5 -0.5 -0.5];
max_mapper = [0.5 -0.5  0.5  0.5];
if num_inputs==2
    [Weights, layers] = Weights_for_2(type);
elseif num_inputs==3
    [Weights, layers] = Weights_for_3(type);
else
    if mod(num_inputs,2)==0
        number = num_inputs/2;
        Weights{1} = sparse(kron(eye(number), mapper));
        layers{1} = cell(4*number, 1);
        layers{1}(:) = {'poslin'};
        if strcmp(type,'min')
            Weights{2} = sparse(kron(eye(number), min_mapper));
        elseif strcmp(type,'max')
            Weights{2} = sparse(kron(eye(number), max_mapper));
        end
        [Weights_remainder , layers_remainder]  =  pure_poslin(type, number);
    else
        number = floor(num_inputs/2);
        Weights{1} = sparse(blkdiag(kron(eye(number), [1 1; -1 -1 ; 1 -1 ; -1 1]) , 1));
        layers{1} = cell(4*number+1, 1);
        layers{1}(1:4*number) = {'poslin'};
        layers{1}(4*number+1) = {'purelin'};
        if strcmp(type,'min')
            Weights{2} = sparse(blkdiag(kron(eye(number), min_mapper),1));
        elseif strcmp(type,'max')
            Weights{2} = sparse(blkdiag(kron(eye(number), max_mapper),1));
        end
        [Weights_remainder , layers_remainder]  =  pure_poslin(type, number+1);
    end
    for i = 1:length(Weights_remainder)
        if i==1
            Weights{2} = sparse(Weights_remainder{1}*Weights{2});
        else
            Weights{i+1} = sparse(Weights_remainder{i});
        end
        if i< length(Weights_remainder)
            layers{i+1} = layers_remainder{i};
        end
    end
end
end




function [Weights, layers] = Weights_for_2(type)

mapper = [1 1; -1 -1 ; 1 -1 ; -1 1];
min_mapper = [0.5 -0.5 -0.5 -0.5];
max_mapper = [0.5 -0.5  0.5  0.5];
Weights{1} = mapper;
layers{1} = cell(4, 1);
layers{1}(:) = {'poslin'};
if strcmp(type, 'min')
    Weights{2} = min_mapper;
elseif strcmp(type,'max')
    Weights{2} = max_mapper;
end
end



function [Weights, layers] = Weights_for_3(type)

mapper = [1 1; -1 -1 ; 1 -1 ; -1 1];
min_mapper = [0.5 -0.5 -0.5 -0.5];
max_mapper = [0.5 -0.5  0.5  0.5];
Weights{1} = blkdiag([1 1; -1 -1 ; 1 -1 ; -1 1] , 1);
layers{1} = cell(4+1, 1);
layers{1}(1:4) = {'poslin'};
layers{1}(4+1) = {'purelin'};
if strcmp(type, 'min')
    Weights{2} = mapper*blkdiag(min_mapper,1);
elseif strcmp(type,'max')
    Weights{2} = mapper*blkdiag(max_mapper,1);
end
layers{2} = cell(4,1);
layers{2}(:) = {'poslin'};
if strcmp(type, 'min')
    Weights{3} = min_mapper;
elseif strcmp(type,'max')
    Weights{3} = max_mapper;
end
end