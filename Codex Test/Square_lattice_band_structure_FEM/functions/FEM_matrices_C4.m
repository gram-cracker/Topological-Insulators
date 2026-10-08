function [M,Mx1,My1,Mx2,My2,Me,Mu,idx] = FEM_matrices_C4(p,t,t_num,xy_box,eps_inc,mu_inc,eps_bk,mu_bk,pol)
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
%   eps: permittivity of the inclusion
%
%   mu_inv: inverse of permeability tensor of the inclusion
%
%
% Outputs:
%
%   M: double, (sparse), finite element mass matrix int_D a*u'*v' dx
%
%   Mx1: double, (sparse), finite element mass matrix int_D a*(u')_x * v dx
%
%   My1: double, (sparse), finite element mass matrix int_D a*(u')_y * v dx
%
%   Mx2: double, (sparse), finite element mass matrix int_D a*(v')_x * u dx
%
%   My2: double, (sparse), finite element mass matrix int_D a*(v')_y * u dx
%
%   Me: double (sparse), finite element stiffness matrix int_D c * u * v dx
%
%   Mu: double, (sparse), finite element mass matrix int_D a * u*v dx
%
% ------------------------------------------------------------------------------

n_nodes = size(p,1);
[n_elements, n_dof] = size(t);

x_min = xy_box(1); x_max = xy_box(2); 
y_min = xy_box(3); y_max = xy_box(4);

if (pol == 'TE')
    coeff_a = 1/mu_inc; coeff_c = eps_inc;
    coeff_a0 = 1/mu_bk; coeff_c0 = eps_bk;
elseif (pol == 'TM')
    coeff_a = 1/eps_inc;  coeff_c = mu_inc;
    coeff_a0 = 1/eps_bk;  coeff_c0 = mu_bk;
else
    disp('Wrong EM ploarization');
    return;
end


% ==============================================================================
% Label nodes
% ==============================================================================
 tol = 1e-5;   % tolerance should be smaller than mesh width
 idx = zeros(n_nodes,1);
 n_equations = 0;
 
 % label the nodes intially, not including nodes on the right and top boundaries
 
 for i = 1:n_nodes      
     if ( abs(p(i,1)-x_max)<tol & abs(p(i,2)-y_max)<tol )
         % corner node: 0 
         idx(i) = 0;
     elseif ( abs(p(i,1)-x_max)<tol )
         % right bnd nodes: -1
         idx(i) = -1;
     elseif ( abs(p(i,2)-y_max)<tol )
         % top bnd nodes: -2
         idx(i) = -2;
     else
         % interior nodes
         n_equations = n_equations + 1;
         idx(i) = n_equations;
     end
 end


% periodic boundary mapping: map the nodes on right and top boundaries to 
% the nodes on the left and bottom boundaries

corner_node = find( abs(p(:,1)-x_min)<tol & abs(p(:,2)-y_min)<tol ); % left bottom corner
left_node = find(abs(p(:,1)-x_min)<tol & abs(p(:,2)-y_max)>tol);     % left bnd
bot_node = find(abs(p(:,2)-y_min)<tol & abs(p(:,1)-x_max)>tol);      % bottom bnd

idx_left = idx(left_node); idx_bot = idx(bot_node);

for i = 1:n_nodes
    if idx(i) == 0
        idx(i) = idx(corner_node);
    elseif idx(i) == -1
        idx(i) = idx_left( abs(p(i,2)-p(left_node,2))<tol );
    elseif idx(i) == -2
        idx(i) = idx_bot( abs(p(i,1)-p(bot_node,1))<tol );
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

mat  = zeros(n_nonzeros,1);    % ( a*dphi_i, dphi_j )
matx1 = zeros(n_nonzeros,1);   % ( phi_i, (a*dphi_j)_x )
maty1 = zeros(n_nonzeros,1);   % ( phi_i, (a*dphi_j)_x )
matx2 = zeros(n_nonzeros,1);   % ( (mu_inv*dphi_j)_x, phi_i )
maty2 = zeros(n_nonzeros,1);   % ( (mu_inv*dphi_j)_x, phi_i )
mate = zeros(n_nonzeros,1);    % ( eps*phi_i, phi_j )
matu = zeros(n_nonzeros,1);    % ( mu_inv(1,1)*phi_i, phi_j )
ir = zeros(n_nonzeros,1);      % row indices
ic = zeros(n_nonzeros,1);      % column indices

[r,s,w] = twod_gauss(7);       % Gauss nodes and weights on standard element
count = 0;



for ne = 1:n_elements
    
    nodes_loc = t(ne,:);
    xy_loc = p(nodes_loc,:);
    xy_center = sum(xy_loc,1)/3;
    
    % Compute local gauss points/weights and shape functions there
    [~,wg,phi,px,py] = twod_shape(xy_loc,r,s,w);  
    
    % permittivity and permeability
    flag_inc = (t_num(ne) == 1);
    coeff_a_loc = flag_inc * coeff_a + (~flag_inc) * coeff_a0;
    coeff_c_loc = flag_inc * coeff_c + (~flag_inc) * coeff_c0;
    
    % coeff_a_loc = 1;
    % coeff_c_loc = fun_medium_coeff(xy_center(1),xy_center(2),coeff_c,coeff_c0,0.2,0.2);
    
    % if (coeff_c_loc == coeff_c)
    % hold on; plot(xy_center(1),xy_center(2),'ro');
    % elseif (coeff_c_loc == coeff_c0)
    % hold on; plot(xy_center(1),xy_center(2),'go');
    % end
    
    % local stiffness matrix    
    mat_loc = twod_bilinear(coeff_a_loc,px,px,wg) + twod_bilinear(coeff_a_loc,py,py,wg);
    matx1_loc = twod_bilinear(coeff_a_loc,px,phi,wg);
    maty1_loc = twod_bilinear(coeff_a_loc,py,phi,wg);
    matx2_loc = twod_bilinear(coeff_a_loc,phi,px,wg);
    maty2_loc = twod_bilinear(coeff_a_loc,phi,py,wg);
    mate_loc = twod_bilinear(coeff_c_loc,phi,phi,wg);
    matu_loc = twod_bilinear(coeff_a_loc,phi,phi,wg);
    
    for ni = 1:n_dof
        for nj = 1:n_dof
            
            count = count + 1;
            ir(count) = idx(nodes_loc(ni));
            ic(count) = idx(nodes_loc(nj));
            
            mat(count) = mat_loc(ni,nj);
            matx1(count) = matx1_loc(ni,nj);
            maty1(count) = maty1_loc(ni,nj);
            matx2(count) = matx2_loc(ni,nj);
            maty2(count) = maty2_loc(ni,nj);
            mate(count) = mate_loc(ni,nj);
            matu(count) = matu_loc(ni,nj);
             
        end
    end
    
    
end


% keep only nonzero entries
ir = ir(1:count);
ic = ic(1:count);
mat = mat(1:count);
matx1 = matx1(1:count);
maty1 = maty1(1:count);
mate = mate(1:count);
matu = matu(1:count);

M = sparse(ir,ic,mat,n_equations,n_equations);
Mx1 = sparse(ir,ic,matx1,n_equations,n_equations);
My1 = sparse(ir,ic,maty1,n_equations,n_equations);
Mx2 = sparse(ir,ic,matx2,n_equations,n_equations);
My2 = sparse(ir,ic,maty2,n_equations,n_equations);
Me = sparse(ir,ic,mate,n_equations,n_equations);
Mu = sparse(ir,ic,matu,n_equations,n_equations);


   