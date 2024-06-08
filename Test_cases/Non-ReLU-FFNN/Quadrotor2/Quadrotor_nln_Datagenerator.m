function [theInput, theOutput, maxmin] = Quadrotor_nln_Datagenerator(x0_vec, y0_vec, z0_vec,  nn, timestep, normalization) 

n=6;
m=3;

x0 = linspace(x0_vec(1),x0_vec(2), x0_vec(3));

y0 = linspace(y0_vec(1),y0_vec(2), y0_vec(3));

z0 = linspace(z0_vec(1),z0_vec(2), z0_vec(3));


T=55;

v         = [1 1 1 ];
vLim      = [  x0_vec(3)  y0_vec(3)  z0_vec(3)  ]; 

T1        =  x0_vec(3)* y0_vec(3)* z0_vec(3)*T  ;

dim = 3+0;


ready = false;
Nend = n;
N0 = Nend+m;
Input = zeros(N0,T1);
Output=zeros(Nend,T1);
j=0;
while ~ready
    initial=[x0(v(1))  y0(v(2))   z0(v(3))  0   0   0];
    for i=1:T
        j=j+1;
        u=cntr(nn(i), initial');
        [~,in_out] =  ode45(@(t,x)system_cont(t,x,u),[0 timestep],initial);
        in_out=in_out(end,:)';
        Input(:,j) = [initial';u];
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
    maxout= maxin(1:n,1);
    minout= minin(1:n,1);
    theInput = (b-a) * diag(1./ (maxin-minin) ) * ( Input - minin )  + a ;
    theOutput= (b-a) * diag(1./(maxout-minout)) * (Output - minout)  + a ;
elseif normalization==0
    theInput=Input;
    theOutput=Output;
    maxmin='no normalization';
end
end




function y = cntr(net, x )
   
    len=length(net.weights)-1;
    for i=1:len
        x=tanh(net.weights{i}*x+net.biases{i});
    end
    y=net.weights{end}*x+net.biases{end};
    
end
    