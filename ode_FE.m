% Yevonnael Andrew 
% Assignment 3 - Computational Methods in Mechanics
% LUT University

function [sol, time] = ode_FE(f,U_0,dt,T)
    N_t = floor(T/dt);
    u = zeros(N_t+1, length(U_0));
    t = linspace(0, N_t*dt, length(u));
    % below is already modified to process vector ODE
    u(1,:) = U_0; 
    t(1) = 0;
    for n = 1:N_t
        u(n+1,:) = u(n,:) + dt*f(u(n,:), t(n));
    end
    sol = u;
    time = t;
end