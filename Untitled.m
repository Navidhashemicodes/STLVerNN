clc
clear
close all

disp('%%%% In order to run the toolbox you need to install the following toolbox on your system,')
disp('%%%%     - YALMIP  + MOSEK solver , (only required for sampling based verification)')


go = input('Did you install Yalmip and Mosek? (If yes insert 1)   ');



if go==1
    
    addpath(genpath( 'Tree_class'))
    addpath(genpath( 'src'))
    
    disp('%%% we present verification results for three probelem. ')
    disp('%%%    1- simple tanh-activation model.')
    disp('%%%    2- Linear time varying plant.')
    disp('%%%    3- Neural Network controlled quadrotor system.')
    problem_number = input('Please indicate which problem do you want to check. (Please insert the number of the problem) : ');
    clc
    switch problem_number
        
        case 1
            
            
        case 2
            
            
        case 3
            
    end
            
    
end