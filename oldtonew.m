function [L_new,G_new,H_new,l_new,g_new,h_new] = oldtonew(L,G,H,l,g,h,omega_as)
    % Raman Singh
    % New variables in terms of old variables
    L_new = L + (omega_as*(L^3)*H);
    G_new = G;
    H_new = H;
    l_new = l - (omega_as*3*(L^2)*H*l);
    g_new = g;
    h_new = h - (omega_as*(L^3)*l);
end