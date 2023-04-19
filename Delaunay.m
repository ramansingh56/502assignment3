%% Advanced Orbital Mechanics HW 3 Problem 2
% Raman Singh
% Propogating the equations of motion in terms of Delaunay variables

close all; clear; clc;

% gravitational parameters
mu = 1;

% initial orbital elements
a = 1;
e = 0.5;
i = pi/4;

% time
t0 = 0;
tf = 100;
tspan = linspace(t0,tf,5000);

% initial conditions
[L0,G0,H0] = del(mu,a,e,i);
l0 = 0; g0 = 0; h0 = 0;
x0 = [L0,G0,H0,l0,g0,h0];

% integration
opts_ode = odeset('RelTol',1e-12,'AbsTol',1e-12); % ode
[t,pq] = ode45(@(t,pq) hamil(t,pq), tspan, x0, opts_ode);

for i = 1:length(pq)
    
    % orbital elements from Delaunay variables
    a = pq(i,1)^2;
    e = sqrt(1-((pq(i,2))/(pq(i,1)))^2);
    inc = acos((pq(i,3))/(pq(i,2)));

    % cartesian positions from orbital elements
    [~,r_vec,~,~] = oe2rv(a,e,inc,pq(i,6),pq(i,5),pq(i,4));
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

function dpqdt = hamil(t,pq)
    
    % angular speed
    omega_as = 0.01;

    % differentiation of the action L wrt time
    dpqdt(1) = 0;

    % differentiation of the action G wrt time
    dpqdt(2) = 0;

    % differentiation of the action H wrt time
    dpqdt(3) = 0;

    % differentiation of the angle l wrt time
    dpqdt(4) = 1/(pq(1)^3);

    % differentiation of the angle g wrt time
    dpqdt(5) = 0;

    % differentiation of the angle h wrt time
    dpqdt(6) = omega_as;
    
    dpqdt = dpqdt';

end