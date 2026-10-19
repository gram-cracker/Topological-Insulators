function [M,Mx,My,Mu,idx] = FEM_matrices_parallelgram_PEC_TE_2PC_edge_state(p,t,cell_v,n_sup,n_inc,xy_inc_left,xy_inc_right)

% Description: compute the mass and stiffness matrics for the joint
% structure, each celll consists of three inclusion.
%
% Inputs:
%
%   p: double, (2,n_nodes) array of mesh nodes
%
%   t: int, (n_elements, n_dof) triangular element matrix
%
%   n_sup: number of the super cells on each side
%
%   n_inc: number of the inclusions in each cell
%
%   xy_inc_left/right: the coordiate of the boundary for all inclusions on the left and right
%
% Outputs:
%
%   M: double, (sparse), finite element mass matrix int_D u'*v' dx
%
%   Mx: double, (sparse), finite element mass matrix int_D (u')_x * v dx
%
%   My: double, (sparse), finite element mass matrix int_D (u')_y * v dx
%
%   Mu: double, (sparse), finite element mass matrix int_D u*v dx
%
% ------------------------------------------------------------------------------

n_nodes = size(p,1);
[n_elements, n_dof] = size(t);

% number of nodes on each inclusion boundary
N_inc0_left = size(xy_inc_left,1)/(n_sup*n_inc);
N_inc0_right = size(xy_inc_right,1)/(n_sup*n_inc);


% ==============================================================================

 % label the inclusion boundary nodes
 idx_inc_bnd = zeros(size(p,1),1);
 for n = 1:n_sup
     for j = 1:n_inc
         n_start_left = (n-1)*(N_inc0_left*n_inc) + (j-1)*N_inc0_left+1;
         n_end_left = (n-1)*(N_inc0_left*n_inc) + j*N_inc0_left;
         [~,label_on_left_inc] = inpolygon(p(:,1),p(:,2),xy_inc_left(n_start_left:n_end_left,1),xy_inc_left(n_start_left:n_end_left,2));
         idx_inc_bnd = idx_inc_bnd + label_on_left_inc;
     end
     for j = 1:n_inc
         n_start_right = (n-1)*(N_inc0_right*n_inc) + (j-1)*N_inc0_left+1;
         n_end_right = (n-1)*(N_inc0_right*n_inc) + j*N_inc0_left;
         [~,label_on_right_inc] = inpolygon(p(:,1),p(:,2),xy_inc_right(n_start_right:n_end_right,1),xy_inc_right(n_start_right:n_end_right,2));
         idx_inc_bnd = idx_inc_bnd + label_on_right_inc;
     end     
 end
 
% Label nodes
% ==============================================================================
 tol = 1e-3;   % tolerance should be smaller than mesh width
 idx = zeros(n_nodes,1);
 n_equations = 0;
 
 % slope of the cell and the lattice constant
 slope = (cell_v(4,2)-cell_v(1,2))/(cell_v(4,1)-cell_v(1,1));
 lat_c = norm(cell_v(4,:)-cell_v(1,:),2);
 
 % the lattic vector and the center of the reference cell
 a1 = [ lat_c, 0 ];
 a2 = cell_v(4,:) - cell_v(1,:);
 % xy_c0 = (a1+a2)/2;

 
 % label the nodes intially, not including nodes on the boundaries on the top/right
 
 for i = 1:n_nodes      
     if ( abs(p(i,2)-slope*(p(i,1)-cell_v(1,1))) < tol | abs(p(i,2)-slope*(p(i,1)-cell_v(2,1))) < tol )
         % left/right bnd nodes: 0; Dirichlet bnd cond
         % hold on; plot(p(i,1),p(i,2),'ro');
         idx(i) = 0;
     elseif ( abs(p(i,2)-cell_v(4,2)) < tol )
         % top bnd nodes: -1; Periodic bnd cond
         % hold on; plot(p(i,1),p(i,2),'go');
         idx(i) = -1;
     else
                  
         if ( idx_inc_bnd(i) > 0 )
             % inclusion bnd nodes
             idx(i) = 0;
             hold on; plot(p(i,1),p(i,2),'-ro')
         else
             % interior nodes
             n_equations = n_equations + 1;
             idx(i) = n_equations;
         end           
         
     end
 end


% periodic boundary mapping: map the nodes on top boundary to 
% the nodes on the bottom boundary

bot_node = find(abs(p(:,2) - cell_v(1,2)) < tol);
% hold on; plot(p(bot_node,1),p(bot_node,2),'go');
idx_bot = idx(bot_node);


% map the nodes on the top boundary to the nodes on the bottom boundary
for i = 1:n_nodes
    if (idx(i) == -1)
        
        temp = p(i,:) - a2 - p(bot_node,:);
        idx(i) = idx_bot( abs(temp(:,1))<tol & abs(temp(:,2))<tol );
        % j = bot_node(abs(temp(:,1))<tol & abs(temp(:,2))<tol);
        % hold on; plot(p(i,1),p(i,2),'go');
        % hold on; plot(p(j,1),p(j,2),'go');
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


mat  = zeros(n_nonzeros,1);    % ( dphi_i, dphi_j )
matx = zeros(n_nonzeros,1);    % ( phi_i, (dphi_j)_x )
maty = zeros(n_nonzeros,1);    % ( phi_i, (dphi_j)_x )
matu = zeros(n_nonzeros,1);    % ( phi_i, phi_j )
ir = zeros(n_nonzeros,1);      % row indices
ic = zeros(n_nonzeros,1);      % column indices

[r,s,w] = twod_gauss(7);       % Gauss nodes and weights on standard element
count = 0;


for ne = 1:n_elements
    
    nodes_loc = t(ne,:);
    xy_loc = p(nodes_loc,:);
    
    % Compute local gauss points/weights and shape functions there
    [~,wg,phi,px,py] = twod_shape(xy_loc,r,s,w);
    
    % local stiff matrix    
    mat_loc = twod_bilinear(1,px,px,wg) + twod_bilinear(1,py,py,wg);
    matx_loc = twod_bilinear(1,px,phi,wg);
    maty_loc = twod_bilinear(1,py,phi,wg);
    matu_loc = twod_bilinear(1,phi,phi,wg);
    
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
matx = matx(1:count);
maty = maty(1:count);
matu = matu(1:count);

M = sparse(ir,ic,mat,n_equations,n_equations);
Mx = sparse(ir,ic,matx,n_equations,n_equations);
My = sparse(ir,ic,maty,n_equations,n_equations);
Mu = sparse(ir,ic,matu,n_equations,n_equations);
   