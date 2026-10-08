function [M,Mx1,My1,Mx2,My2,Me,Mu,idx] = FEM_matrices_C4_2PC_vert(p,t,t_num,xy_box,eps_inc_top,eps_inc_bot,eps_bk,mu_inc_top,mu_inc_bot,mu_bk,pol)
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
%   M: double, (sparse), finite element mass matrix int_D a * u'*v' dx
%
%   Mx1: double, (sparse), finite element mass matrix int_D a * (u')_x * v dx
%
%   My1: double, (sparse), finite element mass matrix int_D a * (u')_y * v dx
%
%   Mx2: double, (sparse), finite element mass matrix int_D u * a * (v')_x  dx
%
%   My2: double, (sparse), finite element mass matrix int_D u * a * (u')_y  dx
%
%   Me: double (sparse), finite element stiffness matrix int_D c * u'*v' dx
%
%   Mu: double, (sparse), finite element mass matrix int_D u*v dx
%
% ------------------------------------------------------------------------------

n_nodes = size(p,1);
[n_elements, n_dof] = size(t);

x_min = xy_box(1); x_max = xy_box(2); 
y_min = xy_box(3); y_max = xy_box(4);

if (pol == 'TE')
    coeff_a_top = 1/mu_inc_top; coeff_c_top = eps_inc_top;
    coeff_a_bot = 1/mu_inc_bot; coeff_c_bot = eps_inc_bot;
    coeff_a0 = 1/mu_bk;         coeff_c0 = eps_bk;
elseif (pol == 'TM')
    coeff_a_top = 1/eps_inc_top;  coeff_c_top = mu_inc_top;
    coeff_a_bot = 1/eps_inc_bot;  coeff_c_top = mu_inc_bot;
    coeff_a0 = 1/eps_bk;          coeff_c0 = mu_bk;
else
    display('Wrong EM ploarization');
    return;
end

% ==============================================================================
% Label nodes
% ==============================================================================
 tol = 1e-3;   % tolerance should be smaller than mesh width
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



% Assembly
% ==============================================================================

n_nonzeros = n_dof^2*n_elements;  % estimate number of non-zero entries

mat  = zeros(n_nonzeros,1);    % ( a*dphi_i, dphi_j )
matx1 = zeros(n_nonzeros,1);   % ( phi_i, (a*dphi_j)_x )
maty1 = zeros(n_nonzeros,1);   % ( phi_i, (a*dphi_j)_y )
matx2 = zeros(n_nonzeros,1);   % ( (a*dphi_i)_x, phi_j,  )
maty2 = zeros(n_nonzeros,1);   % ( (a*dphi_i)_y, phi_j )
mate = zeros(n_nonzeros,1);    % ( c*phi_i, phi_j )
matu = zeros(n_nonzeros,1);    % ( a*phi_i, phi_j )
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
    flag_inc_top = (t_num(ne) == 1);
    flag_inc_bot = (t_num(ne) == 2);
    flag_vac = (t_num(ne) == 3);
    
    % permittivity and permeability on the left and right
    coeff_a_loc = flag_inc_top * coeff_a_top + flag_inc_bot * coeff_a_bot + flag_vac * coeff_a0;
    coeff_c_loc = flag_inc_top * coeff_c_top + flag_inc_bot * coeff_c_bot + flag_vac * coeff_c0;
    
      
    if ( abs(coeff_c_loc - 1) < 1e-3 )
        hold on; plot(xy_center(1),xy_center(2),'ro');
    elseif (abs(coeff_c_loc - 11.7) < 1e-3 )
        hold on; plot(xy_center(1),xy_center(2),'go');
    elseif (abs(coeff_c_loc - 50) < 1e-3 )
        hold on; plot(xy_center(1),xy_center(2),'co');
    end
    
    % local stiff matrix    
    mat_loc = twod_bilinear(coeff_a_loc,px,px,wg) + twod_bilinear(coeff_a_loc,py,py,wg);
    matx1_loc = twod_bilinear(coeff_a_loc,px,phi,wg);
    maty1_loc = twod_bilinear(coeff_a_loc,py,phi,wg);
    matx2_loc = twod_bilinear(coeff_a_loc,phi,px,wg);
    maty2_loc = twod_bilinear(coeff_a_loc,phi,py,wg);
    mate_loc = twod_bilinear(coeff_c_loc,phi,phi,wg);
    matu_loc = twod_bilinear(coeff_a_loc,phi,phi,wg);
    
    for ni = 1:n_dof
        if (idx(nodes_loc(ni))>0) 
            for nj = 1:n_dof
                if (idx(nodes_loc(nj))>0)
                    
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
    end
    
    
end


% keep only nonzero entries
ir = ir(1:count);
ic = ic(1:count);
mat = mat(1:count);
matx1 = matx1(1:count);
maty1 = maty1(1:count);
matx2 = matx2(1:count);
maty2 = maty2(1:count);
mate = mate(1:count);
matu = matu(1:count);

M = sparse(ir,ic,mat,n_equations,n_equations);
Mx1 = sparse(ir,ic,matx1,n_equations,n_equations);
My1 = sparse(ir,ic,maty1,n_equations,n_equations);
Mx2 = sparse(ir,ic,matx2,n_equations,n_equations);
My2 = sparse(ir,ic,maty2,n_equations,n_equations);
Me = sparse(ir,ic,mate,n_equations,n_equations);
Mu = sparse(ir,ic,matu,n_equations,n_equations);


   