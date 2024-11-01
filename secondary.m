syms F_f phi r_s_min d_s k_s kappa_s d0s theta0s T0 T1 tau_s T1s1 T1s2 T0s1 T0s2 r_helix theta_helix N_helix F_bs

% theta_s = d_s/(r_helix*tan(theta_helix));
% tau_ss = kappa_s*(theta0s + theta_s);
% F_ss = k_s*(d0s + d_s);
r_s = d_s/tan(phi) + r_s_min;

syms theta_s tau_ss F_ss

eq1 = 0 == -tau_s - tau_ss + T1s1*r_s - T0s1*r_s + r_helix*N_helix*sin(theta_helix);
eq2 = 0 == tau_ss + T1s2*r_s - T0s2*r_s - r_helix*N_helix*sin(theta_helix);
eq3 = 0 == F_ss - F_bs + N_helix*cos(theta_helix);
eq4 = T1 == T1s1 + T1s2;
eq5 = T0 == T0s1 + T0s2;
eq6 = F_f == T1 - T0;
eq7 = T1s1 == T1s2;
eq8 = T0s1 == T0s2;

S = solve([eq1 eq2 eq3 eq4 eq5 eq6 eq7 eq8], [d_s T1s1 T0s1 T1s2 T0s2 T1 N_helix F_bs])

S_ds = simplify(S.d_s)
S_N_helix = simplify(S.N_helix)
S_F_bs = simplify(S.F_bs)

% 0 == F_ss - F_bs + tau_ss/r_helix + F_helix
% F_helix = 