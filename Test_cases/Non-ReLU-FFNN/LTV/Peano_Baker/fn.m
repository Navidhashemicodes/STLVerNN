function dydt = fn(t,y)
    dydt = zeros(2,1);
    u = zeros(2,1);
    u(1) = sin(y(1)) - cos(y(2));
    u(2) = 0;

    dydt(1) = y(2) + u(1);
    dydt(2) = -y(2) - (2 + sin(t))*y(1) + u(2);
end

% u(1) = sin(y(1)) + y(2)^2
% u(2) = y(1) + y(2)
