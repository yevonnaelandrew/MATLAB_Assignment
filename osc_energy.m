% Yevonnael Andrew 
% Assignment 3 - Computational Methods in Mechanics
% LUT University

% This function calculates the kinetic energy, the potential energy, 
% and return the total energy as the sum of KE and PE.

function TE = osc_energy(u,v,omega)
    PE = 0.5*omega^2*u^2;
    KE = 0.5*v^2;
    TE = PE+KE;
end