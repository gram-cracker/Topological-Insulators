function [ LDOS, grad_LDOS ] = fun_LDOS(epsilon,kappa,omega,mesh,FEM_matrices,src,F)
% my_fun_LDOS computes the LDOS functions and the gradient


% polarization of EM wave
pol = 'TE';

% desired frequency omega_0 and frequency band gamma
omega0 = real(omega);
gamma = imag(omega);

% source locations and strength
xy_src = src.xy;
sign_src = src.sign;


% -------------------------------------------------------------------------
% the mesh data
% -------------------------------------------------------------------------
p = mesh.p;
t = mesh.t;
t_num = mesh.t_num;
idx = mesh.idx;
n_elements = size(mesh.t,1);
% area = zeros(n_elements,1);

% ---------------------------------------------------------------------
% FEM matrics
% ---------------------------------------------------------------------
Me = FEM_matrices.Me;
M = FEM_matrices.M;
Mx1 = FEM_matrices.Mx1;
My1 = FEM_matrices.My1;
Mx2 = FEM_matrices.Mx2;
My2 = FEM_matrices.My2;
Mu = FEM_matrices.Mu;

A = M - 1i * ( kappa(1)*Mx1 + kappa(2)*My1 ) +  1i * ( kappa(1)*Mx2 + kappa(2)*My2 ) + (kappa*kappa')*Mu;
A_conj = M + 1i * ( kappa(1)*Mx1 + kappa(2)*My1 ) -  1i * ( kappa(1)*Mx2 + kappa(2)*My2 ) + (kappa*kappa')*Mu;
B = Me;
Opt = A - omega^2 * B;
Opt_conj = A_conj - omega^2 * B;

% ---------------------------------------------------------------------
% Forward problem
% ---------------------------------------------------------------------
u = Opt\(1i*omega*F);  u_conj = Opt_conj\(1i*omega*F);
uu = u(idx);  uu_conj = u_conj(idx);

% ---------------------------------------------------------------------
% LDOS
% ---------------------------------------------------------------------   
u_src = u_interp( uu, p, t, xy_src(:,1), xy_src(:,2) );
uconj_src = u_interp( uu_conj, p, t, xy_src(:,1), xy_src(:,2) );
u_dot_uconj = inner_prod3(p,t,t_num,epsilon,pol,uu,uu_conj);
LDOS = real(-6/pi * (sign_src * u_src) + 6/pi * (1i*gamma)/(omega0+1i*gamma) * sign_src * uconj_src ) ...
         + real(12*gamma/pi * u_dot_uconj);
% LDOS = real(-6/pi * (sign_src * u_src) + 6/pi * (1i*gamma)/(omega0+1i*gamma) * sign_src * uconj_src);
% LDOS = 12*gamma/pi * real(u_dot_uconj);

% ---------------------------------------------------------------------
% The gradient of LDOS
% --------------------------------------------------------------------- 

if (nargout > 1)
    % Adjoint problem
    src_adj1 = - Me * u_conj;
    u_adj1 = Opt_conj\src_adj1;
    uu_adj1 = u_adj1(idx);

    src_adj2 = - Me * u;
    u_adj2 = Opt\src_adj2;
    uu_adj2 = u_adj2(idx);

    % compute the gradient
    grad_LDOS = zeros(n_elements,1);
    for ne = 1:n_elements
        nodes_loc = t(ne,:);
        
        xi = p(nodes_loc(1),1); yi = p(nodes_loc(1),2); 
        xj = p(nodes_loc(2),1); yj = p(nodes_loc(2),2); 
        xk = p(nodes_loc(3),1); yk = p(nodes_loc(3),2); 
        area_e = det([1 xi yi; 1 xj yj; 1 xk yk])/2;
        
        % area(ne) = area_e;

        uu_loc = sum(uu(nodes_loc))/3;
        uu_conj_loc = sum(uu_conj(nodes_loc))/3;
        uu_adj1_loc = sum(uu_adj1(nodes_loc))/3;
        uu_adj2_loc = sum(uu_adj2(nodes_loc))/3;
        grad_LDOS(ne) = area_e * real( 6/pi * (1i*omega0+2*gamma) * (uu_loc .* uu_conj_loc)  ...
                        - 12*gamma*omega^2/pi * (uu_loc .* uu_adj1_loc + uu_conj_loc .* uu_adj2_loc) );                                 
        % grad_LDOS(ne) = area_e * real(6/pi * (1i*omega0) * (uu_loc .* uu_conj_loc));
        % grad_LDOS(ne) = area_e * real(12*gamma/pi * (uu_loc * uu_conj_loc) - 12*gamma*omega^2/pi * (uu_loc * uu_adj1_loc + uu_conj_loc * uu_adj2_loc));

    end
    
end


end
