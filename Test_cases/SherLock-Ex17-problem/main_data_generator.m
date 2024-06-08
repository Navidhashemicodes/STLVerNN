clear all
clc
close all
x_vec=[ 0.35-0.3   0.4+0.3 30];
y_vec=[-0.35-0.3  -0.3+0.3 30];
z_vec=[ 0.35-0.3   0.4+0.3 30];

nn_controller =  NN_Reader( 3,1,'neural_network_controller');

normalization=1;
timestep=0.2;

[theInput, theOutput, maxmin] = Ex17_nln_Datagenerator(x_vec, y_vec, z_vec, nn_controller, timestep, normalization);

Input=theInput;
Output=theOutput;
save('Data.mat','Input', 'Output');
maxin=maxmin.maxin;
minin=maxmin.minin;
maxout=maxmin.maxin(1:3,1);
minout=maxmin.minin(1:3,1);
save('maxmin.mat','maxin', 'maxout', 'minin', 'minout');

clear all
