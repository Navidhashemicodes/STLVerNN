function [theInput, theOutput, maxmin] = ACC_nln_Datagenerator(pos1_vec, pos2_vec, veloc1_vec, veloc2_vec, Hidst1_vec, Hidst2_vec, nn, timestep, normalization) 

n=6;
m=1;

pos1=linspace(pos1_vec(1),pos1_vec(2), pos1_vec(3));

pos2=linspace(pos2_vec(1),pos2_vec(2), pos2_vec(3));

veloc1=linspace(veloc1_vec(1),veloc1_vec(2), veloc1_vec(3));

veloc2=linspace(veloc2_vec(1),veloc2_vec(2), veloc2_vec(3));

Hidst1=linspace(Hidst1_vec(1),Hidst1_vec(2), Hidst1_vec(3));

Hidst2=linspace(Hidst2_vec(1),Hidst2_vec(2), Hidst2_vec(3));

E=[ 0  0  0  0  1  0;...
    1  0  0 -1  0  0;...
    0  1  0  0 -1  0];

V_set=30;
t_gap=1.4;
T=50;

v         = [1 1 1 1 1 1 ];
vLim      = [pos1_vec(3)  veloc1_vec(3) Hidst1_vec(3)  pos2_vec(3)  veloc2_vec(3)  Hidst2_vec(3) ];   
T1        =  pos1_vec(3)* pos2_vec(3)* veloc1_vec(3)* veloc2_vec(3)* Hidst1_vec(3)* Hidst2_vec(3)*T  ;
dim=6+0;


ready = false;
Nend=n;
N0=Nend+m;
Input=zeros(N0,T1);
Output=zeros(Nend,T1);
j=0;
while ~ready
    initial=[pos1(v(1))  veloc1(v(2))  Hidst1(v(3)) pos2(v(4)) veloc2(v(5)) Hidst2(v(6))];
    for i=1:T
        j=j+1;
        init_a=[V_set; t_gap; E*initial'];
        a_ego=cntr(nn, init_a);
        [~,in_out] =  ode45(@(t,x)dynamicsACC(t,x,a_ego),[0 timestep],initial);
        in_out=in_out(end,:)';
        Input(:,j) = [initial';a_ego];
        Output(:,j)= in_out;
        initial= in_out';
    end
    
   
    ready = true;
    ff=v(dim);
    for k = 1:dim                                   %index updater
        v(k) = v(k) + 1;
        if v(k) <= vLim(k)
            ready = false;
            break;
        end
        v(k) = 1;
    end
    
   
end

if normalization==1
    a=-1;
    b=1;
    maxin = max(Input')';
    maxmin.maxin=maxin;
    minin = min(Input')';
    maxmin.minin=minin;
    maxout= maxin(1:6,1);
    minout= minin(1:6.1);
    theInput = (b-a) * diag(1./ (maxin-minin) ) * ( Input - minin )  + a ;
    theOutput= (b-a) * diag(1./(maxout-minout)) * (Output - minout)  + a ;
elseif normalization==0
    theInput=Input;
    theOutput=Output;
    maxmin='no normalization';
end
end




function y = cntr(net, x)
   
    len=length(net.weights)-1;
    for i=1:len
        x=poslin(net.weights{i}*x+net.biases{i});
    end
    y=net.weights{end}*x+net.biases{end};
    
end
    