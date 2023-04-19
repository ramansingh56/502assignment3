function [L,G,H] = del(mu,a,e,i)

    % mean motion
    n = sqrt(mu/a^3);

    % delaunay variables
    L = n*a^2;
    G = L*sqrt(1-e^2);
    H = G*cos(i);

end