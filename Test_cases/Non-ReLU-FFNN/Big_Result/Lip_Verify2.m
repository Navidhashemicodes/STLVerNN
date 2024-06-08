function [decision, Letsgo] = Lip_Verify2( net, border, lb,ub, M )
x_axis=linspace(lb(1),ub(1),M+1);
y_axis=linspace(lb(2),ub(2),M+1);
Center=cell(M,M);
epsilon=cell(M,M);
eps= norm(ub-lb,2)/2;
decision = 'Yes';
Letsgo=1;
for i=1:M
    for j=1:M
        clear alpha_params beta_params  Lb  Ub
        Center{i,j}= 0.5*([x_axis(i);y_axis(j)]+ [x_axis(i+1);y_axis(j+1)]);
        epsilon{i,j}= 0.5*([x_axis(i+1);y_axis(j+1)]- [x_axis(i);y_axis(j)]);
        rho2(1)=NN(net, Center{i,j}-epsilon{i,j});
        rho2(2)=NN(net, Center{i,j}+epsilon{i,j});
        rho2(3)=NN(net, Center{i,j}+[1;-1].*epsilon{i,j});
        rho2(4)=NN(net, Center{i,j}+[-1;1].*epsilon{i,j});
        Rho2(i,j)=min(rho2);
        if Rho2(i,j)<0
            decision= 'No';
            Letsgo=0;
        end
        method_Trajectory='Linear-CROWN';
%         method_Trajectory='Quadratic-CROWN';
%         method_Trajectory='Quadratic-Gurobi';
        method_STL       ='approx-star' ;
%         method_STL       ='exact-star' ;
        [Lb, Ub, alpha_params, beta_params]= Trapezius_bounds_General(Center{i,j}, epsilon{i,j}, net, border,  method_Trajectory, method_STL);
        LipSDP_type= 'LipSDP-layer';
        addpath C:\Program Files\Mosek
        [Lip(i,j), status(i,j), time(i,j)] = Trapezius_Lip_SDP( net, border, alpha_params, beta_params , Lb, Ub,  LipSDP_type );
        
        certificate= Rho2(i,j)/Lip(i,j);
        
        if eps>certificate
            [decision, Letsgo] = Lip_Verify(net,border, Center{i,j}-epsilon{i,j}, Center{i,j}+epsilon{i,j}, M);
        end
        save1=Lip(i,j);
        save2=Rho2(i,j);
        save3=Center{i,j};
        save4=epsilon{i,j};
        save([num2str(save3(1)) '_' num2str(save3(2)) '_' num2str(save4(1)) '.mat']   , 'save1' , 'save2' , 'save3' , 'save4')
        if Letsgo==0
            break;
        end
    end
    if Letsgo==0
        break;
    end
end