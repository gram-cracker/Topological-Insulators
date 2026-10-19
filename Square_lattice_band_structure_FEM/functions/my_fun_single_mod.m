function [ f, df ] = my_fun_single_mod(epsilon,kappa,omega,xy_box,mesh,FEM_matrices,src,F)

% -------------------------------------------------------------------------
% the mesh data
% -------------------------------------------------------------------------
p = mesh.p;
t = mesh.t;
t_num = mesh.t_num;

mu = ones(size(t,1),1);
Me = mass_matrix_C4(p,t,t_num,xy_box,epsilon,mu,'TE');
FEM_matrices.Me = Me;

[ LDOS, grad_LDOS ] = fun_LDOS(epsilon,kappa,omega,mesh,FEM_matrices,src,F);

f = -LDOS;
df = -grad_LDOS;

end

