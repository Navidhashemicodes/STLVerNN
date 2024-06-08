clear all
clc
close all
pos1_vec=[80 110 120];
pos2_vec=[10 11 10];
veloc1_vec=[32 32.2 5];
veloc2_vec=[30 30.2 5];
Hidst1_vec=[0 0 1];
Hidst2_vec=[0 0 1];
load('controller_3_20.mat')
nn.weights=network.W;
nn.biases =network.b;
normalization=1;
timestep=0.1;

[theInput, theOutput, maxmin] = ACC_nln_Datagenerator(pos1_vec, pos2_vec, veloc1_vec, veloc2_vec, Hidst1_vec, Hidst2_vec, nn, timestep, normalization);

Input=theInput;
Output=theOutput;
save('Data.mat','Input', 'Output', 'maxmin');

clear all
