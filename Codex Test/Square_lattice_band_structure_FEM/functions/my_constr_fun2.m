function [c,ceq,DC,DCeq] = my_constr_fun2(epsilon,kappa,omega,xy_box,mesh,FEM_matrices,src_even,src_odd,F_even,F_odd)

dkappa = 0.1;
slope = 0.1;
n_eig = [ 4,5];

[f,g] = my_fun_multi_mod_multi_kappa(epsilon(1:end-1),kappa,omega,xy_box,mesh,FEM_matrices,src_even,src_odd,F_even,F_odd);
[f_slope_p,g_slope_p] = my_fun_omega_slope(epsilon(1:end-1),kappa,dkappa,n_eig,xy_box,mesh,FEM_matrices);
[f_slope_m,g_slope_m] = my_fun_omega_slope(epsilon(1:end-1),kappa,-dkappa,n_eig,xy_box,mesh,FEM_matrices);

f_slope_p
f_slope_m

% c = f - epsilon(end);
c = [ f - epsilon(end), f_slope_p(1)+slope, slope-f_slope_p(2), slope-f_slope_m(1), f_slope_m(2)+slope ];
ceq = [];
% DC = [g;-1,-1,-1,-1,-1,-1];
DC = [g, g_slope_p(:,1), -g_slope_p(:,2), -g_slope_m(:,1), g_slope_m(:,2); -1,-1,-1,-1,-1,-1,0,0,0,0];
DCeq = [];

end