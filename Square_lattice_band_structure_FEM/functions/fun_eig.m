function [ omega, grad_omega ] = fun_eig(epsilon,kappa,n_eig,mesh,FEM_matrices)
% my_fun_eig computes n_eig'th eigenvalues and the gradients w.r.t epsilon


max_num_eig = max(n_eig);

% -------------------------------------------------------------------------
% the mesh data
% -------------------------------------------------------------------------
p = mesh.p;
t = mesh.t;
t_num = mesh.t_num;
idx = mesh.idx;
n_elements = size(mesh.t,1);

% ---------------------------------------------------------------------
% FEM matrics
% ---------------------------------------------------------------------
M = FEM_matrices.M;
Mx1 = FEM_matrices.Mx1;
My1 = FEM_matrices.My1;
Mx2 = FEM_matrices.Mx2;
My2 = FEM_matrices.My2;
Mu = FEM_matrices.Mu;
Me = FEM_matrices.Me;


A = M - 1i * ( kappa(1)*Mx1 + kappa(2)*My1 ) +  1i * ( kappa(1)*Mx2 + kappa(2)*My2 ) + (kappa*kappa')*Mu;
B = Me;



% ---------------------------------------------------------------------
% eigenvalue and the corresponding eigenfunction
% ---------------------------------------------------------------------   
% solve the eigenvalues
[ U, lambda ] = eigs(A,B,max_num_eig,'smallestabs','Tolerance',1e-5);
omega_all = sqrt(diag(lambda));
omega = omega_all(n_eig)';


if (nargout > 1)  
    
    grad_omega = zeros(n_elements,length(n_eig));

    for m = 1:length(n_eig)
    
        % the eigenmode
        u = U(:,n_eig(m));
        uu = u(idx);
        norm_uu2 = inner_prod3(p,t,t_num,epsilon,'TE',uu,conj(uu));
        uu = uu/sqrt(norm_uu2);

        % ---------------------------------------------------------------------
        % The gradient of omega
        % ---------------------------------------------------------------------          
        for ne = 1:n_elements
        
            nodes_loc = t(ne,:);
        
            xi = p(nodes_loc(1),1); yi = p(nodes_loc(1),2); 
            xj = p(nodes_loc(2),1); yj = p(nodes_loc(2),2); 
            xk = p(nodes_loc(3),1); yk = p(nodes_loc(3),2); 
            area_e = det([1 xi yi; 1 xj yj; 1 xk yk])/2;
      
            % uu_loc = sum(uu(nodes_loc))/3;
            % grad_omega(ne) = -area_e * omega/2 * (abs(uu_loc))^2;
        
            uu_loc = uu(nodes_loc);
            mat_loc = area_e/12 * [ 2 1 1; 1 2 1; 1 1 2 ];
            grad_omega(ne,m) = -omega(m)/2 * (uu_loc'*mat_loc*uu_loc);
            
        end
    
    end
    
end


omega = real(omega);
grad_omega = real(grad_omega);

end
