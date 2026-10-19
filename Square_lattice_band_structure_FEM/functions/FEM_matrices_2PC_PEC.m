function [M,Mx,My,Me,Mu,idx] = FEM_matrices_2PC_PEC(p,t,t_num,xy_box,xc_top,xc_bot,R)
% Description: compute the mass and stiffness matrics
%
% Inputs:
%
%   p: double, (2,n_nodes) array of mesh nodes
%
%   t: int, (n_elements, n_dof) triangular element matrix
%
%   t_num: array of regions numbers for each triangle in t
%
%   xy_box: double, [x_min, x_max, y_min, y_max] denote the domain
%
%   xc_top, xc_bot: centers of the inclusions in the top/bottom cells
%
%   R: radius of inclusions

% Outputs:
%
%   M: double, (sparse), finite element mass matrix int_D mu_inv*u'*v' dx
%
%   Mx: double, (sparse), finite element mass matrix int_D (mu_inv*u')_x * v dx
%
%   My: double, (sparse), finite element mass matrix int_D (mu_inv*u')_y * v dx
%
%   Me: double (sparse), finite element stiffness matrix int_D e u'*v' dx
%
%   Mu: double, (sparse), finite element mass matrix int_D mu_inv(1,1) * u*v dx
%
% ------------------------------------------------------------------------------

n_nodes = size(p,1);
[n_elements, n_dof] = size(t);

x_min = xy_box(1); x_max = xy_box(2); 
y_min = xy_box(3); y_max = xy_box(4);

% ==============================================================================
% Label nodes
% ==============================================================================
 tol = 1e-5;   % tolerance should be smaller than mesh width
 idx = zeros(n_nodes,1);
 n_equations = 0;
 
 % label the nodes intially, not including nodes on the right and top boundaries
 
 for i = 1:n_nodes
     if ( ( abs(p(i,2)-y_min)<tol ) | ( abs(p(i,2)-y_max)<tol ) )
         % top or bottom bnd nodes: 0 (Dirichelt bnd cond)
         idx(i) = 0;
     elseif ( abs(p(i,1)-x_max)<tol )
         % right bnd nodes: -1
         idx(i) = -1;
     elseif ( min(sqrt( (p(i,1)-xc_top(:,1)).^2 + (p(i,2)-xc_top(:,2)).^2 )) <= R+1e-10 ...
              | min(sqrt( (p(i,1)-xc_bot(:,1)).^2 + (p(i,2)-xc_bot(:,2)).^2 )) <= R+1e-10   )
                   % inclusion bnd nodes
         idx(i) = 0;
         hold on; plot(p(i,1),p(i,2),'-o')
     else
         % interior nodes
         n_equations = n_equations + 1;
         idx(i) = n_equations;
     end
 end

% periodic boundary mapping: map the nodes on right boundary to 
% the nodes on the left boundary

left_node = find(abs(p(:,1)-x_min)<tol & abs(p(:,2)-y_min)>tol & abs(p(:,2)-y_max)>tol);    
idx_left = idx(left_node);

for i = 1:n_nodes
    if idx(i) == -1        
        idx(i) = idx_left( abs(p(i,2)-p(left_node,2))<tol );
    end
end

if ( ~isempty(find(idx<0)) )
    disp('periodic boundary mapping not successful, adjust mesh size');
    return
end


% ==============================================================================
% Assembly
% ==============================================================================

n_nonzeros = n_dof^2*n_elements;  % estimate number of non-zero entries

mat  = zeros(n_nonzeros,1);    % ( mu*dphi_i, dphi_j )
matx = zeros(n_nonzeros,1);    % ( phi_i, (mu*dphi_j)_x )
maty = zeros(n_nonzeros,1);    % ( phi_i, (mu*dphi_j)_x )
mate = zeros(n_nonzeros,1);    % ( eps*phi_i, phi_j )
matu = zeros(n_nonzeros,1);    % ( mu_inv(1,1)*phi_i, phi_j )
ir = zeros(n_nonzeros,1);      % row indices
ic = zeros(n_nonzeros,1);      % column indices

[r,s,w] = twod_gauss(7);       % Gauss nodes and weights on standard element
count = 0;


for ne = 1:n_elements
    
    nodes_loc = t(ne,:);
    xy_loc = p(nodes_loc,:);
    % xy_center = sum(xy_loc,1)/3;
    
    % Compute local gauss points/weights and shape functions there
    [~,wg,phi,px,py] = twod_shape(xy_loc,r,s,w);  
    
    % local stiff matrix    
    mat_loc = twod_bilinear(1,px,px,wg) + twod_bilinear(1,py,py,wg);
    matx_loc = twod_bilinear(1,px,phi,wg);
    maty_loc = twod_bilinear(1,py,phi,wg);
    mate_loc = twod_bilinear(1,phi,phi,wg);
    
    for ni = 1:n_dof
        if (idx(nodes_loc(ni))>0) 
            for nj = 1:n_dof
                if (idx(nodes_loc(nj))>0)
            
                    count = count + 1;
                    ir(count) = idx(nodes_loc(ni));
                    ic(count) = idx(nodes_loc(nj));
            
                    mat(count) = mat_loc(ni,nj);
                    matx(count) = matx_loc(ni,nj);
                    maty(count) = maty_loc(ni,nj);
                    mate(count) = mate_loc(ni,nj);
                end
            
            end
        end
    end
    
    
end


% keep only nonzero entries
ir = ir(1:count);
ic = ic(1:count);
mat = mat(1:count);
matx = matx(1:count);
maty = maty(1:count);
mate = mate(1:count);
matu = matu(1:count);

M = sparse(ir,ic,mat,n_equations,n_equations);
Mx = sparse(ir,ic,matx,n_equations,n_equations);
My = sparse(ir,ic,maty,n_equations,n_equations);
Me = sparse(ir,ic,mate,n_equations,n_equations);
Mu = sparse(ir,ic,matu,n_equations,n_equations);


   