clc
clear
close all

disp('%%%% In order to run the toolbox you need to install the following toolbox on your system,')
disp('%%%%     - NNV toolbox: you can download it by  command "  git clone https://github.com/transafeailab/nnv.git  " ')


go = input('Did you install NNV? (If yes insert 1)   ');



if go==1
    
    addpath(genpath( 'Tree_class'))
    addpath(genpath( 'src'))
    
    tabel_generator = input('Please indicate which table do you want to check. (Please insert the number of the table) : ');
    clc
    switch tabel_generator
        
        case 1
            
            problem_number = input('Please indicate the number of STL property. (As an example insert 1 for property $\varphi_1$):  ');
            clc
            switch problem_number
                case 1
                    disp('%%%  This property is verified on a ReLU-NN model of dimension [3,10,10,10,2] using a controller of')
                    disp('%%%  dimension [2,50,1,2,1,2,1,1]. This property is verified with reachability technique.')
                    disp('%%% We utilize Approx-star and Exact-star for the verification.' )
                    Reach_method = input('Please indicate which technique do you want to use for verification. (For Exact-star insert 1 and For Approx-star insert 2)  :  ');
                    switch Reach_method
                        case 1
                            cd Test_cases\SherLock-Ex14-problem
                            Test_phi1_exact
                            cd ..
                            cd ..
                            clc
                            disp(['The Run-time of verification algorithm is ' num2str(Computation_time)]);
                            disp(['The Robustness interval from the verification algorithm is [' num2str(interval)  ']. ']);
                            disp( 'The verification was succesfull since lower-bound is positive.');
                        case 2
                            cd Test_cases\SherLock-Ex14-problem
                            Test_phi1_Approx
                            cd ..
                            cd ..
                            clc
                            disp(['The Run-time of verification algorithm is ' num2str(Computation_time)]);
                            disp(['The Robustness interval from the verification algorithm is [' num2str(interval)  ']. ']);
                            disp( 'Unfortunately the verification failed since lower-bound is negative.');
                            
                    end
                case 2
                    disp('%%%  This property is verified on a ReLU-NN model of dimension [4,10,10,3] using a ReLU NN controller')
                    disp('%%%  of dimension [3,100,1,2,1,2,1,1]. This property is verified with reachability technique.')
                    disp('%%%  We utilize Approx-star and Exact-star for the verification.' )
                    Reach_method = input('Please indicate which technique do you want to use for verification. (For Exact-star insert 1 and For Approx-star insert 2)  :  ');
                    switch Reach_method
                        case 1
                            cd Test_cases\SherLock-Ex17-problem
                            Test_phi2_exact
                            cd ..
                            cd ..
                            clc
                            disp(['The Run-time of verification algorithm is ' num2str(Computation_time)]);
                            disp(['The Robustness interval from the verification algorithm is [' num2str(interval)  ']. ']);
                            disp( 'The verification was succesfull since lower-bound is positive.');
                        case 2
                            cd Test_cases\SherLock-Ex17-problem
                            Test_phi2_Approx
                            cd ..
                            cd ..
                            clc
                            disp(['The Run-time of verification algorithm is ' num2str(Computation_time)]);
                            disp(['The Robustness interval from the verification algorithm is [' num2str(interval)  ']. ']);
                            disp( 'Unfortunately the verification failed since lower-bound is negative.');
                            
                    end
                    
                case 3
                    disp('%%%  This property is verified in two cases on Adaptive cruise control problem.')
                    disp('%%%  The first case assumes a linear model and the second addresses a nonlinear model with friction.')
                    nln = input( 'In order to verify the property for the non-linear model press 1. Otherwise press 2 for linear model.');
                    clc
                    
                    switch nln
                        case 1
                            disp('%%%  This property is verified on a ReLU-NN model of dimension [7,10,10,6] using a ReLU NN controller')
                            disp('%%%  of dimension [5,20,20,20,1]. This property is verified with reachability technique.')
                            disp('%%%  We utilize Approx-star and Exact-star for the verification.' )
                            Reach_method = input('Please indicate which technique do you want to use for verification. (For Exact-star insert 1 and For Approx-star insert 2)  :  ');
                            clc
                            switch Reach_method
                                case 1
                                    cd Test_cases\SherLock-ACC-problem
                                    Test_phi3_exact
                                    cd ..
                                    cd ..
                                    clc
                                    disp(['The Run-time of verification algorithm is ' num2str(Computation_time)]);
                                    disp(['The Robustness interval from the verification algorithm is [' num2str(interval)  ']. ']);
                                    disp( 'The verification was succesfull since lower-bound is positive.');
                                case 2
                                    cd Test_cases\SherLock-ACC-problem
                                    Test_phi3_Approx
                                    cd ..
                                    cd ..
                                    clc
                                    disp(['The Run-time of verification algorithm is ' num2str(Computation_time)]);
                                    disp(['The Robustness interval from the verification algorithm is [' num2str(interval)  ']. ']);
                                    disp( 'The verification was succesfull since lower-bound is positive.');
                                    
                            end
                            
                        case 2
                            disp('%%%  This property is verified on a LTI model of dimension [7,6] using a ReLU NN controller')
                            disp('%%%  of dimension [5,20,20,20,1]. This property is verified with reachability technique.')
                            disp('%%%  We utilize Approx-star and Exact-star for the verification.' )
                            Reach_method = input('Please indicate which technique do you want to use for verification. (For Exact-star insert 1 and For Approx-star insert 2)  :  ');
                            clc
                            switch Reach_method
                                case 1
                                    cd Test_cases\Trapezius_for_ACC
                                    Test_phi3_linear_exact
                                    cd ..
                                    cd ..
                                    clc
                                    disp(['The Run-time of verification algorithm is ' num2str(Computation_time)]);
                                    disp(['The Robustness interval from the verification algorithm is [' num2str(interval)  ']. ']);
                                    disp( 'The verification was succesfull since lower-bound is positive.');
                                case 2
                                    cd Test_cases\Trapezius_for_ACC
                                    Test_phi3_linear_Approx
                                    cd ..
                                    cd ..
                                    clc
                                    disp(['The Run-time of verification algorithm is ' num2str(Computation_time)]);
                                    disp(['The Robustness interval from the verification algorithm is [' num2str(interval)  ']. ']);
                                    disp( 'The verification was succesfull since lower-bound is positive.');
                                    
                            end
                    end
                    
                case 4
                    disp('%%%  This property is verified on a ReLU-NN model of dimension [4,8,2] using a ReLU NN controller')
                    disp('%%%  of dimension [2,8,2]. This property is verified with reachability technique.')
                    disp('%%%  We utilize Exact-star for the verification.' )
                    verify = input('please press 1 and then hit enter to start the verification:  ');
                    switch verify
                        case 1
                            cd Test_cases\STL_Complexity
                            Test_until1_new
                            cd ..
                            cd ..
                            clc
                            disp(['The Run-time of verification algorithm is ' num2str(Computation_time)]);
                            disp(['The Robustness interval from the verification algorithm is [' num2str(interval)  ']. ']);
                            disp( 'The verification was succesfull since lower-bound is positive.');
                        case 2
                            
                    end
                    
                case 5
                    disp('%%%  This property is verified on a ReLU-NN model of dimension [4,8,2] using a ReLU NN controller')
                    disp('%%%  of dimension [2,8,2]. This property is verified with reachability technique.')
                    disp('%%%  We utilize Exact-star for the verification.' )
                    verify = input('please press 1 and the hit enter to start the verification:  ');
                    switch verify
                        case 1
                            cd Test_cases\STL_Complexity
                            Test_until2_new
                            cd ..
                            cd ..
                            clc
                            disp(['The Run-time of verification algorithm is ' num2str(Computation_time)]);
                            disp(['The Robustness interval from the verification algorithm is [' num2str(interval)  ']. ']);
                            disp( 'The verification was succesfull since lower-bound is positive.');
                        case 2
                            
                    end
                case 6
                    disp('%%%  This property is verified on a ReLU-NN model of dimension [4,8,2] using a ReLU NN controller')
                    disp('%%%  of dimension [2,8,2]. This property is verified with reachability technique.')
                    disp('%%%  We utilize Exact-star for the verification.' )
                    verify = input('please press 1 and then hit enter to start the verification:  ');
                    switch verify
                        case 1
                            cd Test_cases\STL_Complexity
                            Test_until3_new
                            cd ..
                            cd ..
                            clc
                            disp(['The Run-time of verification algorithm is ' num2str(Computation_time)]);
                            disp(['The Robustness interval from the verification algorithm is [' num2str(interval)  ']. ']);
                            disp( 'The verification was succesfull since lower-bound is positive.');
                        case 2
                            
                    end
                case 7
                    disp('%%%  This property is successfully rejected on a ReLU-NN model of dimension [4,8,2] using a ReLU NN controller')
                    disp('%%%  of dimension [2,8,2]. This property is rejected with reachability technique.')
                    disp('%%%  We utilize Exact-star for the verification.' )
                    verify = input('please press 1 and then hit enter to start the verification:  ');
                    switch verify
                        case 1
                            cd Test_cases\STL_Complexity
                            Test_until4_new
                            cd ..
                            cd ..
                            clc
                            disp(['The Run-time of verification algorithm is ' num2str(Computation_time)]);
                            disp(['The Robustness interval from the verification algorithm is [' num2str(interval)  ']. ']);
                            disp( 'The Rejection of property was succesfull since lower-bound is negative and the algorithm is sound and complete.');
                        case 2
                            
                    end
            end
            
            
            
            
        case 2
            disp('%%% In this table we verify the exponential stability for an LTI model with dimension [2,2], the model')
            disp('%%% is controlled with a ReLU NN controller of dimension [2, 30, 30, 30, 2]. The property is verified')
            disp('%%% with reachability technique.')
            disp('%%% We utilize exact-star for the verification.')
            disp('%%% We partitioned the set of initial condition in 25 different sub-region. ')
            clear all
            clc
            close all
            center = input('Please specify the center of the sub-region which you are interested in and is presented on the table (For example: [-41;88] ):  ');
            verify = input('please press 1 and then hit enter to start the verification:  ');
            switch verify
                case 1
                    cd Test_cases\Exponential_Stability
                    [Computation_time , interval ] = exponential_stability( center);
                    cd ..
                    cd ..
                    clc
                    disp(['The Run-time of verification algorithm is ' num2str(Computation_time)]);
                    disp(['The Robustness interval from the verification algorithm is [' num2str(interval)  ']. ']);
                    disp( 'The verification of property was succesfull since lower-bound is positive. ');
                case 2
                    
            end
            
    end
    
    
else
    disp('The codes will not run if NNV is not installed')
end

    
