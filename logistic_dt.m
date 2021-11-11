% Yevonnael Andrew 
% Assignment 3 - Computational Methods in Mechanics
% LUT University

% This function help to find the best timesteps for a particular logistic
% function. At each iteration, it will plot the last two timesteps and 
% ask user whether to continue to halve or to end the iteration.

f = @(u,t) 0.1*(1 - u/500)*u; % logistic function
U_0 = 100; % intial value
dt = 50; % timesteps
T = 200; % periods
k = 1; % helper iterator

% the function will iterate indefinitely until user terminate the loop
while true
    % using Forward Euler with the last two timesteps
    [u1,t1] = ode_FE(f, U_0, dt/(2^k), T);
    [u2,t2] = ode_FE(f, U_0, dt/(2^(k-1)), T);
    plot(t1,u1,'b-',t2,u2,'r--'); % plot the last two timesteps
    xlabel(['k: ', num2str(k), ' blue (dt_k): ', num2str(dt/(2^k)), ', red (dt_k_-_1): ', num2str(dt/(2^(k-1)))]); ylabel('N(t)');
    % asking user whether to continue or to terminate
    is_cont = input("Continue? 'y'/'n': ", 's');
    if is_cont == 'y'
        k = k+1;
    else 
        break;
    end
end