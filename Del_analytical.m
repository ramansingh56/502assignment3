%% Advanced Orbital Mechanics HW 3 Problem 2
% Raman Singh
% Analytical - Delaunay variables

close all; clear; clc;

% time
t0 = 0;
tf = 100;
tspan = linspace(t0,tf,5000);

for i = 1:length(tspan)
    [L,G,H,l,g,h] = ana_del(tspan(i));
    
    % orbital elements from Delaunay variables
    a = L^2;
    e = sqrt(1-(G/L)^2);
    inc = acos(H/G);

    % cartesian positions from orbital elements
    [~,r_vec,~,~] = oe2rv(a,e,inc,h,g,l);
    rad_vec(i,:) = r_vec;
end

% plotting
plot3(rad_vec(:,1),rad_vec(:,2),rad_vec(:,3))
hold on
plot3(0,0,0,'r*')
xlabel('x in DU')
ylabel('y in DU')
zlabel('z in DU')
title('Trajectory plot')

function [L,G,H,l,g,h] = ana_del(t)
    
    % gravitational parameters
    mu = 1;
    
    % initial orbital elements
    a = 1;
    e = 0.5;
    i = pi/4;
    
    % initial conditions
    [L0,G0,H0] = del(mu,a,e,i);
    l0 = 0; g0 = 0; h0 = 0;

    % angular speed
    omega_as = 0.01;

    L = L0;

    G = G0;

    H = H0;
    
    l = (t/(L^3)) +l0;

    g = g0;

    h = omega_as*t+h0;

end