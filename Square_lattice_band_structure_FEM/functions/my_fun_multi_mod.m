function [ f, df ] = my_fun_multi_mod(epsilon,kappa,omega,xy_box,mesh,FEM_matrices,src1,src2,F1,F2)                                      

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
[ LDOS1, grad_LDOS1 ] = fun_LDOS(epsilon,kappa,omega,mesh,FEM_matrices,src1,F1);
[ LDOS2, grad_LDOS2 ] = fun_LDOS(epsilon,kappa,omega,mesh,FEM_matrices,src2,F2);

f = -[LDOS1, LDOS2];
df = -[ grad_LDOS1, grad_LDOS2];

end

