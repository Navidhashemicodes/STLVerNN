clear
clc
close all

xxs = -3:0.01:3;
yys = -3:0.01:3;

normalization=0;

index=0;
for xx = 1:length(xxs)
    for yy = 1:length(yys)
        index=index+1;
        s = [xxs(xx) yys(yy)];
        Input(:,index) = s;
        Output(:,index) = sin(s(1))-cos(s(2));
    end
end


if normalization==1
    a=-1;
    b=1;
    maxin = max(Input')';
    maxmin.maxin=maxin;
    minin = min(Input')';
    maxmin.minin=minin;
    maxout= max(Output')';
    maxmin.maxout=maxout;
    minout= min(Output')';
    maxmin.minout=minout;
    theInput = (b-a) * diag(1./ (maxin-minin) ) * ( Input - minin )  + a ;
    theOutput= (b-a) * diag(1./(maxout-minout)) * (Output - minout)  + a ;
elseif normalization==0
    theInput=Input;
    theOutput=Output;
    maxmin='no normalization';
end



Input=theInput;
Output=theOutput;
save('Data.mat','Input', 'Output');

if normalization==1
    maxin=maxmin.maxin;
    minin=maxmin.minin;
    maxout=maxmin.maxout;
    minout=maxmin.minout;
    save('maxmin.mat','maxin', 'maxout', 'minin', 'minout');
end
clear all