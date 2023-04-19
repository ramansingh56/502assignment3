%% Advanced Orbital Mechanics Problem 1 & 2
% Raman Singh
% Comparative study between the fully perturbed solution and the 2nd order
% accurate perturbed solution

close all; clear; clc;

% gravitational parameter
mu = 1;

% initial orbital elements
a_0 = 0.1:0.1:1.5;
e_0 = 0.1:0.1:1;
i_0 = pi/18:pi/18:4*pi/9;

% angular speed
omega_as = 0.5;

a_c = randsample(a_0,1);
e_c = randsample(e_0,1);
i_c = randsample(i_0,1);

% time
t0 = 0;
tf = 100;
tspan = linspace(t0,tf,5000);

% initial conditions for old variables
[L0,G0,H0] = del(mu,a_c,e_c,i_c);
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
    e = sqrt(1-((pq(i,2))/(pq(i,1)))^2);
    inc = acos((pq(i,3))/(pq(i,2)));
    
    % equinoctial elements
    he_ham(i) = e*sin(g+h);
    ke_ham(i) = e*cos(g+h);
    pe_ham(i) = tan(inc/2)*sin(h);
    qe_ham(i) = tan(inc/2)*cos(h);

end

%% Kamiltonian approach - 2nd order approximate solution

% integration
[tk,pqn] = ode45(@(tk,pqn) kamil(tk,pqn), tspan, x_new0, opts_ode);

for i = 1:length(pqn)
    
    [L,G,H,l,g,h] = newtoold(pqn(i,1),pqn(i,2),pqn(i,3),pqn(i,4),pqn(i,5),pqn(i,6),omega_as);
    e = sqrt(1-((G/L)^2));
    inc = acos(H/G);

    % equinoctial elements
    he_kam(i) = e*sin(g+h);
    ke_kam(i) = e*cos(g+h);
    pe_kam(i) = tan(inc/2)*sin(h);
    qe_kam(i) = tan(inc/2)*cos(h);

end

% plotting
figure
subplot(2,1,1)
plot(he_ham,ke_ham,Color='k',LineWidth=1.5)
hold on
plot(he_kam,ke_kam,Color='r',LineWidth=1)
xlabel('h')
ylabel('k')
title('h vs. k')

subplot(2,1,2)
plot(pe_ham,qe_ham,Color='k',LineWidth=1.5)
hold on
plot(pe_kam,qe_kam,Color='r',LineWidth=1)
xlabel('p')
ylabel('q')
title('p vs. q')
sgtitle("\omega = " + omega_as + " TU^{-1}, a = " + a_c + ", e = " + e_c + ", i = " + i_c + "")

% legend('Fully perturbed','2nd order accurate')

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