clear all
clc
close all

for i=1:64
load(['Results_hh_' num2str(i-1) '.mat']); T(i)=run_time/60;
end


bar(T)
hold on
linee=mean(T)*ones(1,66);
plot(0:65, linee);