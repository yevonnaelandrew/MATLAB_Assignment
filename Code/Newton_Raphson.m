% Yevonnael Andrew 
% Assignment 2 of 2 - Computational Methods in Mechanics
% LUT University

% Initialize time array to loop
t_a=linspace(0,1,101);
% Initialize empty array to store results
theta_array=zeros(101,1);
d_array=zeros(101,1);
theta_dot_array=zeros(101,1);
d_dot_array=zeros(101,1);

% Looping over time array
for t=1:length(t_a)
    % Solve F and J to obtain theta and d
    [x, ~] = Newton_system(@F, @J, [pi/4,0.2], t_a(t), 0.0001);
    theta_array(t)=x(1); % theta
    d_array(t)=x(2); % d
    % Solve derivative of F and J to obtain theta_dot and d_dot
    [x_dot, ~] = Newton_time_derivative(@F_t, @J_t, [-0.4,0.05], t_a(t), theta_array(t), 0.0001);
    theta_dot_array(t)=x_dot(1); % theta_dot
    d_dot_array(t)=x_dot(2); % d_dot
end

% Plot theta, distance, time and their derivative
figure; plot(t_a, theta_array, t_a, d_array);
xlabel('Time')
legend({'theta', 'd'})
title('Time vs Theta & Distance Plot')
figure; plot(t_a, theta_dot_array, t_a, d_dot_array);
xlabel('Time')
legend({'theta', 'd'})
title('Time Derivative Plot')

% Taken from book
function [x, iteration_counter] = Newton_system(F, J, x, t, eps)
    % Solve nonlinear system F=0 by Newton’s method.
    % J is the Jacobian of F. Both F and J must be functions of x.
    % At input, x holds the start value. The iteration continues
    % until ||F|| < eps.
    F_value = F(x,t);
    F_norm = norm(F_value); % l2 norm of vector
    iteration_counter = 0;
    while abs(F_norm) > eps && iteration_counter < 100
        delta = J(x,t)\-F_value;
        x = x + delta;
        F_value = F(x,t);
        F_norm = norm(F_value);
        iteration_counter = iteration_counter + 1;
    end
    % Here, either a solution is found, or too many iterations
    if abs(F_norm) > eps
        iteration_counter = -1;
    end
end

% Modified version of above function to solve the time derivative
function [x, iteration_counter] = Newton_time_derivative(F_t, J_t, x, t, theta, eps)
    % Solve nonlinear system F=0 by Newton’s method.
    % J is the Jacobian of F. Both F and J must be functions of x.
    % At input, x holds the start value. The iteration continues
    % until ||F|| < eps.
    F_value = F_t(x,t,theta);
    F_norm = norm(F_value); % l2 norm of vector
    iteration_counter = 0;
    while abs(F_norm) > eps && iteration_counter < 100
        delta = J_t(t)\-F_value;
        x = x + delta;
        F_value = F_t(x,t,theta);
        F_norm = norm(F_value);
        iteration_counter = iteration_counter + 1;
    end
    % Here, either a solution is found, or too many iterations
    if abs(F_norm) > eps
        iteration_counter = -1;
    end
end

function F_matrix = F(x,t)
    a = 0.1;
    b = 0.2;
    phi = pi/6 - 1*t;
    theta = x(1);
    d = x(2);
    F_matrix = [a*cos(phi) + b*cos(theta) - d; ...
                a*sin(phi) - b*sin(theta)];
end

% Jacobian of F
function J_matrix = J(x,t)
    b = 0.2;
    theta = x(1);
    J_matrix = [-b*sin(theta) -1; ...
                -b*cos(theta) 0];
end

% Derivative of F with respect to time
function F_t_matrix = F_t(x,t,theta)
    a = 0.1;
    b = 0.2;
    phi = pi/6 - 1*t;
    phi_dot = -1;
    theta_dot=x(1);
    d_dot=x(2);
    F_t_matrix = [-a*phi_dot*sin(phi) - b*theta_dot*sin(theta) - d_dot; ...
                   a*phi_dot*cos(phi) - b*theta_dot*cos(theta)];
end

% Jacobian of F_t
function J_t_matrix = J_t(theta)
    b = 0.2;
    J_t_matrix = [-b*sin(theta) -1; ...
                  -b*cos(theta) 0];
end