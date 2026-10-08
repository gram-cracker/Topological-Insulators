function [ f, df ] = my_fun_multi_mod_multi_kappa(epsilon,kappa,omega,xy_box,mesh,FEM_matrices,src1,src2,F1,F2)                                      

kappa_shift = [ 0.01, 0 ];
omega_shift = 0.02;

% -------------------------------------------------------------------------
% the mesh data
% -------------------------------------------------------------------------
p = mesh.p;
t = mesh.t;
t_num = mesh.t_num;

% -------------------------------------------------------------------------
% the FEM matrices
% -------------------------------------------------------------------------
mu = ones(size(t,1),1);
Me = mass_matrix_C4(p,t,t_num,xy_box,epsilon,mu,'TE');
FEM_matrices.Me = Me;

% -------------------------------------------------------------------------
% the LDOS
% -------------------------------------------------------------------------
[ LDOS1_minus, grad_LDOS1_minus ] = fun_LDOS(epsilon,kappa-kappa_shift,omega-omega_shift,mesh,FEM_matrices,src1,F1);
[ LDOS2_minus, grad_LDOS2_minus ] = fun_LDOS(epsilon,kappa-kappa_shift,omega+omega_shift,mesh,FEM_matrices,src2,F2);
[ LDOS1, grad_LDOS1 ] = fun_LDOS(epsilon,kappa,omega,mesh,FEM_matrices,src1,F1);
[ LDOS2, grad_LDOS2 ] = fun_LDOS(epsilon,kappa,omega,mesh,FEM_matrices,src2,F2);
[ LDOS1_plus, grad_LDOS1_plus ] = fun_LDOS(epsilon,kappa+kappa_shift,omega+omega_shift,mesh,FEM_matrices,src1,F1);
[ LDOS2_plus, grad_LDOS2_plus ] = fun_LDOS(epsilon,kappa+kappa_shift,omega-omega_shift,mesh,FEM_matrices,src2,F2);

f = -[LDOS1_minus, LDOS2_minus, LDOS1, LDOS2, LDOS1_plus, LDOS2_plus];
df = -[ grad_LDOS1_minus, grad_LDOS2_minus, grad_LDOS1, grad_LDOS2, grad_LDOS1_plus, grad_LDOS2_plus];

end