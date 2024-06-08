clear all
clc
close all
x_vec=[0.8-0.3 0.9+0.3 100];
y_vec=[0.4-0.3 0.5+0.3 100];

nn_controller =  NN_Reader( 2,1,'neural_network_controller');

normalization=1;
timestep=0.1;

[theInput, theOutput, maxmin] = Ex14_nln_Datagenerator(x_vec, y_vec,  nn_controller, timestep, normalization);

Input=theInput;
Output=theOutput;
save('Data.mat','Input', 'Output');
maxin=maxmin.maxin;
minin=maxmin.minin;
maxout=maxmin.maxin(1:2,1);
minout=maxmin.minin(1:2,1);
save('maxmin.mat','maxin', 'maxout', 'minin', 'minout');

clear all
