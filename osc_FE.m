% Yevonnael Andrew 
% Assignment 3 - Computational Methods in Mechanics
% LUT University

% This function calculates solutions using the Forward-Euler method and
% will return a plot to compare numerical and exact solutions, 
% and a plot of total energy in the oscillation system.

% initial condition
X_0 = 2;
omega = 2;
P = 2*pi/omega;
dt = P/20;
T = 5*P;

[u,v,t,e] = osc_FE_func(X_0, omega, dt, T);
figure; plot(t,u,'b-', t, X_0*cos(omega*t), 'r--');
title('Numerical vs exact');
figure; plot(t,e,'b-');
title('Total Energy with the Forward Euler Method');

function [u,v,t,e] = osc_FE_func(X_0, omega, dt, T)
    N_t=floor(T/dt);
    t=linspace(0,N_t*dt,N_t+1);

    % empty arrays to store values
    u=zeros(N_t+1,1);
    v=zeros(N_t+1,1);
    e=zeros(N_t+1,1);

    u(1) = X_0;
    v(1) = 0;
    e(1) = osc_energy(u(1),v(1),omega);
    
    % iteratively foward with time
    for n = 1:N_t
        % Follow eq 4.47 and 4.48
        u(n+1) = u(n) + dt*v(n);
        v(n+1) = v(n) - dt*omega^2*u(n);
        e(n+1) = osc_energy(u(n+1),v(n+1),omega); % calculate energy
    end
end