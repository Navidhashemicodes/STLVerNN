function [theInput, theOutput, maxmin] = Ex13_nln_Datagenerator(x_vec, y_vec, nn, timestep, normalization) 

n=2;
m=1;

x=linspace(x_vec(1),x_vec(2), x_vec(3));

y=linspace(y_vec(1),y_vec(2), y_vec(3));


T=50;

v         = [1 1];
vLim      = [x_vec(3)  y_vec(3) ];   
T1        =  x_vec(3)* y_vec(3)*T  ;
dim=2+0;


ready = false;
Nend=n;
N0=Nend+m;
Input=zeros(N0,T1);
Output=zeros(Nend,T1);
j=0;
while ~ready
    initial=[x(v(1))  y(v(2)) ];
    for i=1:T
        j=j+1;
        u=cntr(nn, initial');
        [~,in_out] =  ode45(@(t,x)ex_13(t,x,u),[0 timestep],initial);
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
    maxout= maxin(1:2,1);
    minout= minin(1:2.1);
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
    