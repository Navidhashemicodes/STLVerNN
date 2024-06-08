clear all
clc
close all
x0_vec=[  0.0638-0.03   0.1063+0.03 40];
y0_vec=[  0.0638-0.03   0.1063+0.03 40];
z0_vec=[ -0.0213-0.01   0.0213+0.01 42];


n=6;


load('nn_controller.mat')

for i=1:55
    NN_actors=NN_actor;
    weights1 = NN_actor.weights{1};
    biases1 = NN_actor.biases{1};
    biases1 = biases1+(i-1)*weights1(:,7);
    weights1 = weights1(:,1:6);
    NN_actors.weights{1} = weights1;
    NN_actors.biases{1} = biases1;
    nn(i).weights = NN_actors.weights;
    nn(i).biases = NN_actors.biases;
    nn(i).layers = NN_actor.layers;
end


normalization=1;
timestep=0.05;

[theInput, theOutput, maxmin] = Quadrotor_nln_Datagenerator(x0_vec, y0_vec, z0_vec, nn, timestep, normalization);

Input=theInput;
Output=theOutput;
save('Data.mat','Input', 'Output');
if normalization==1
    maxin=maxmin.maxin;
    minin=maxmin.minin;
    maxout=maxmin.maxin(1:n,1);
    minout=maxmin.minin(1:n,1);
    save('maxmin.mat','maxin', 'maxout', 'minin', 'minout');
end

clear all
