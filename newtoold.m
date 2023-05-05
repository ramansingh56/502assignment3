function [L,G,H,l,g,h] = newtoold(L_new,G_new,H_new,l_new,g_new,h_new,omega_as)
    % Raman Singh
    % Old variables in terms of new variables
    L = L_new + (omega_as*((L_new)^3)*H_new);
    G = G_new;
    H = H_new;
    l = l_new - (omega_as*3*((L_new)^2)*H_new*l_new);
    g = g_new;
    h = h_new - (omega_as*((L_new)^3)*l_new);
end