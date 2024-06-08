xxs =-1+2*rand(1,10);
yys = -1+2*rand(1,10);

figure;
hold on;
for xx = 1:length(xxs)
    for yy = 1:length(yys)
        state = [xxs(xx) yys(yy)];
        [tout,yout] = ode45(@fn,[0 10],state);
        plot(tout,yout(:,1),'r',tout,yout(:,2),'b');
        hold on; 
    end
end
