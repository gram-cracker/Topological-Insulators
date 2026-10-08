function [f,g] = my_fun_omega_slope(epsilon,kappa,dkappa,n_eig,xy_box,mesh,FEM_matrices)
%MY_FUN_OMEGA_SLOPE returns the funcions and their gradients on the slopes
% of the eigenfrequency

% perturbation of kappa
% delta = 0.1;
dir = [ 1, 0 ];
kappa_pert = kappa + dkappa * dir;

% mesh and mass matrices
p = mesh.p;
t = mesh.t;
t_num = mesh.t_num;
mu = ones(size(t,1),1);
Me = mass_matrix_C4(p,t,t_num,xy_box,epsilon,mu,'TE');
FEM_matrices.Me = Me;

% compute the eigvalues and their gradients w.r.t. epsilon, at kappa and kappa_pert
[ omega, grad_omega ] = fun_eig(epsilon,kappa,n_eig,mesh,FEM_matrices);
[ omega_pert, grad_omega_pert ] = fun_eig(epsilon,kappa_pert,n_eig,mesh,FEM_matrices);

f = (omega_pert - omega)/dkappa;
g = (grad_omega_pert - grad_omega)/dkappa;

end

