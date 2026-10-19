function Me = mass_matrix_C4(p,t,t_num,xy_box,epsilon,mu,pol)
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
%   epsilon: permittivity of the medium
%
%   mu: permeability of the medium
%
%
% Output:
%
%   Me: double (sparse), finite element stiffness matrix int_D c * u'*v' dx
%
% ------------------------------------------------------------------------------

n_nodes = size(p,1);
[n_elements, n_dof] = size(t);

x_min = xy_box(1); x_max = xy_box(2); 
y_min = xy_box(3); y_max = xy_box(4);

if (pol == 'TE')
    coeff_a = 1./mu; coeff_c = epsilon;
elseif (pol == 'TM')
    coeff_a = 1./epsilon;  coeff_c = mu;
else
    display('Wrong EM ploarization');
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

mate = zeros(n_nonzeros,1);    % ( eps*phi_i, phi_j )
ir = zeros(n_nonzeros,1);      % row indices
ic = zeros(n_nonzeros,1);      % column indices

[r,s,w] = twod_gauss(7);       % Gauss nodes and weights on standard element
count = 0;



for ne = 1:n_elements
    
    nodes_loc = t(ne,:);
    xy_loc = p(nodes_loc,:);
    
    % Compute local gauss points/weights and shape functions there
    [~,wg,phi,px,py] = twod_shape(xy_loc,r,s,w);
    
    % permittivity and permeability
    coeff_c_loc = coeff_c(ne);       
    
    % local stiffness matrix    
    mate_loc = twod_bilinear(coeff_c_loc,phi,phi,wg);
    
    for ni = 1:n_dof
        for nj = 1:n_dof
            
            count = count + 1;
            ir(count) = idx(nodes_loc(ni));
            ic(count) = idx(nodes_loc(nj));
            
            mate(count) = mate_loc(ni,nj);
             
        end
    end
        
end


% keep only nonzero entries
ir = ir(1:count);
ic = ic(1:count);
mate = mate(1:count);

Me = sparse(ir,ic,mate,n_equations,n_equations);   
