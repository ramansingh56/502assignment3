%% Advanced Orbital Mechanics Problem 1 & 2
% Raman Singh
% Check for consistency

close all; clear; clc;

% gravitational parameter
mu = 1;

% initial orbital elements
a = 1;
e = 0.5;
i = pi/4;

% angular speed
omega_as = 0.1;

% time
t0 = 0;
tf = 100;
tspan = linspace(t0,tf,5000);

% initial conditions for old variables
[L0,G0,H0] = del(mu,a,e,i);
l0 = 0; g0 = 0; h0 = 0;
x0 = [L0,G0,H0,l0,g0,h0];

% initial conditions for new variables
[L_new0,G_new0,H_new0,l_new0,g_new0,h_new0] = oldtonew(L0,G0,H0,l0,g0,h0,omega_as);
x_new0 = [L_new0,G_new0,H_new0,l_new0,g_new0,h_new0];

%% Kamiltonian approach - 1st order accurate solution

% integration
opts_ode = odeset('RelTol',1e-12,'AbsTol',1e-12); % ode
[tk,pqn] = ode45(@(tk,pqn) kamil(tk,pqn), tspan, x_new0, opts_ode);

for i = 1:length(pqn)
    L = pqn(i,1); G = pqn(i,2); H = pqn(i,3);
    l = pqn(i,4); g = pqn(i,5); h = pqn(i,6);
    
    a = L^2;
    e = sqrt(1-(G/L)^2);
    inc = acos(H/G);

    % cartesian positions from orbital elements
    [~,r_vec,~,~] = oe2rv(a,e,inc,h,g,l);
    rad_vec_k(i,:) = r_vec;

end

% plotting
figure
plot3(rad_vec_k(:,1),rad_vec_k(:,2),rad_vec_k(:,3),Color='r',LineWidth=1)
xlabel('x in DU')
ylabel('y in DU')
zlabel('z in DU')
title(['Kamiltonian trajectory plot for \omega = ' num2str(omega_as) ' TU^{-1}'])

function dpqndt = kamil(tk,pqn)
    
    % differentiation of the action L' wrt time
    dpqndt(1) = 0;

    % differentiation of the action G' wrt time
    dpqndt(2) = 0;

    % differentiation of the action H' wrt time
    dpqndt(3) = 0;

    % differentiation of the angle l' wrt time
    dpqndt(4) = 1/(pqn(1)^3);

    % differentiation of the angle g' wrt time
    dpqndt(5) = 0;

    % differentiation of the angle h' wrt time
    dpqndt(6) = 0;
    
    dpqndt = dpqndt';

end