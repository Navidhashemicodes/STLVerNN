clear
clc

vlim=[3 3 3 3 3 3];
dim=length(vlim);
events=prod(vlim);
v=ones(1,dim)
for steps=1:events-1
    Stepp=steps;
    for r=2:dim
        remainder=prod(vlim(r:end));
        v_p(r-1) = floor(Stepp/remainder);
        Stepp = mod(Stepp, remainder);
    end
    v_p(dim)=Stepp;
    v=v_p+1
end