clc
clear
close all


figure_generator = input('Please indicate which figure do you want to be plotted. (Please insert the number of the figure. We do not plot figure 1,8,9) : ');
clc
switch figure_generator
    
    case 2
        cd Test_cases\SherLock-Ex14-problem
        reach_set
        cd ..
        cd ..
        clc
    case 3
        cd Test_cases\SherLock-Ex17-problem
        reach_set
        cd ..
        cd ..
        clc
    case 4
        cd Test_cases\Exponential_Stability
        reachset
        cd ..
        cd ..
        clc
    case 5
        go = input('Did you install Mosek and add it to your MATLAB path? (If yes, please insert 1)   ');
        if go == 1
            disp('%%%% We do not plot the figure, but we compute results in this figure as a tree. We know a recursive algorithm generates a tree. The process takes 12 minutes.')
            disp('%%%% The results will be generated in 4 separate trees')
            disp('%%%%     1- The first tree shows the centers of sub-regions in figure 5')
            disp('%%%%     2- The second tree shows the dimeters of sub-regions in figure 5')
            disp('%%%%     3- The third tree shows the Lipschitz constant of sub-regions in figure 5')
            disp('%%%%     4- The fourth tree shows the minimum of robustness of the vertices in each sub-region in figure 5')
            %             disp('%%%% The files contains the certificates rho1 (named as save1), rho2 (named as save2), center of sub-region (named save3) and diameter of it (named save4)')
            %             disp('%%%% The addressed sub-region is also understandable from the file name.')
            %             disp('%%%% For example the file:  1.1875_1.8125_and_0.0625_0.0625.mat  adresses a sub-region with center [1.1875;1.8125] with diameters [0.0625; 0.0625].')
            verify = input('please press 1 and then hit enter to start the verification:  ');
            if verify ==1
                cd Test_cases\Non-ReLU-FFNN\Main
                Test_LIPSDP
                cd ..
                cd ..
                cd ..
                clc
            end
        end
    case 6
        cd Test_cases\Non-ReLU-FFNN\LTV
        reach_set
        cd ..
        cd ..
        cd ..
        clc
    case 7
        cd Test_cases\Non-ReLU-FFNN\Quadrotor2
        reach_set
        cd ..
        cd ..
        cd ..
        clc
        
    case 10
        cd Test_cases\STL_Complexity
        reach_set_until2
        cd ..
        cd ..
        clc
    case 11
        
        cd Test_cases\Non-ReLU-FFNN\Main
        reach_set
        cd ..
        cd ..
        cd ..
        clc
        
end
        

                        