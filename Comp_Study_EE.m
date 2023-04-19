%% Advanced Orbital Mechanics Problem 1 & 2
% Raman Singh
% Comparative study between the fully perturbed solution and the 2nd order
% accurate perturbed solution

close all; clear; clc;

% gravitational parameter
mu = 1;

% initial orbital elements
a = 1;
e = 0.5;
i = pi/4;

% angular speed
omega_as = 0.01;

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

%% Hamiltonian approach - Fully perturbed solution

% integration
opts_ode = odeset('RelTol',1e-12,'AbsTol',1e-12); % ode
[th,pq] = ode45(@(th,pq) hamil(th,pq,omega_as), tspan, x0, opts_ode);

for i = 1:length(pq)
    
    l = pq(i,4); g = pq(i,5); h = pq(i,6);

    % orbital elements from Delaunay variables
    a = pq(i,1)^2;
    e = sqrt(1-((pq(i,2))/(pq(i,1)))^2);
    inc = acos((pq(i,3))/(pq(i,2)));
    
    % cartesian positions from orbital elements
    [~,r_vec,~,~] = oe2rv(a,e,inc,h,g,l);
    rad_vec_h(i,:) = r_vec;

end

%% Kamiltonian approach - 1st order accurate solution

% integration
[tk,pqn] = ode45(@(tk,pqn) kamil(tk,pqn), tspan, x_new0, opts_ode);

for i = 1:length(pqn)
    
    [L,G,H,l,g,h] = newtoold(pqn(i,1),pqn(i,2),pqn(i,3),pqn(i,4),pqn(i,5),pqn(i,6),omega_as);
    a = L^2;
    e = sqrt(1-(G/L)^2);
    inc = acos(H/G);

    % cartesian positions from orbital elements
    [~,r_vec,~,~] = oe2rv(a,e,inc,h,g,l);
    rad_vec_k(i,:) = r_vec;

end

% plotting
figure
plot3(rad_vec_h(:,1),rad_vec_h(:,2),rad_vec_h(:,3),Color='k',LineWidth=1)
hold on
plot3(rad_vec_k(:,1),rad_vec_k(:,2),rad_vec_k(:,3),Color='r',LineWidth=1)
plot3(0,0,0,'g*')
xlabel('x in DU')
ylabel('y in DU')
zlabel('z in DU')
title(['Trajectory plot for \omega = ' num2str(omega_as) ' TU^{-1}'])
legend('Fully perturbed','2nd order accurate','Center')

function dpqdt = hamil(th,pq,omega_as)

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

function dpqndt = kamil(tk,pqn)
    
    % differentiation of the action L wrt time
    dpqndt(1) = 0;

    % differentiation of the action G wrt time
    dpqndt(2) = 0;

    % differentiation of the action H wrt time
    dpqndt(3) = 0;

    % differentiation of the angle l wrt time
    dpqndt(4) = 1/(pqn(1)^3);

    % differentiation of the angle g wrt time
    dpqndt(5) = 0;

    % differentiation of the angle h wrt time
    dpqndt(6) = 0;
    
    dpqndt = dpqndt';

end