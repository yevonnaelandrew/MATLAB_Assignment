% Yevonnael Andrew 
% Assignment 3 - Computational Methods in Mechanics
% LUT University

% This function calculates solutions using the Backward-Euler method and
% will return a plot to compare numerical and exact solutions, 
% and a plot of total energy in the oscillation system.

% initial condition
X_0 = 2;
omega = 2;
P = 2*pi/omega;
dt = P/200; % try different dt to see the difference
T = 3*P;

[u,v,t,e] = osc_BE_func(X_0, omega, dt, T);
figure; plot(t,u,'b-', t, X_0*cos(omega*t), 'r--');
title('Numerical vs exact');
figure; plot(t,e,'b-');
title('Total Energy with the Backward Euler Method');

% function to solve oscillations using Backward Euler method
function [u,v,t,e] = osc_BE_func(X_0, omega, dt, T)
    N_t=floor(T/dt);
    t=linspace(0,N_t*dt,N_t+1);
    
    % empty arrays to store values
    u=zeros(N_t+1,1);
    v=zeros(N_t+1,1);
    e=zeros(N_t+1,1);
    
    u(1) = X_0;
    v(1) = 0;
    e(1) = osc_energy(u(1),v(1),omega);
    
    % using formula derived in the word document, we then calculate the
    % v(n) and u(n) iteratively
    for n = 2:N_t+1
        v(n) = (v(n-1) - dt*omega^2*u(n-1))/(1+dt^2*omega^2);
        u(n) = u(n-1) + dt*v(n);
        e(n) = osc_energy(u(n),v(n),omega); % calculate energy
    end
end