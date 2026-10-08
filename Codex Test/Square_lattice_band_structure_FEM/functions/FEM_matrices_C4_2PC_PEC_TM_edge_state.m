function [M,Mx,My,Mu,M1,M2,idx] = FEM_matrices_C4_2PC_PEC_TM_edge_state(p,t,xy_box,n_sup,n_inc,xy_c,xy_inc_left,xy_inc_right)
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
%   n_sup: number of the supercells on each side
%
%   n_inc: number of inclusions in each cell
%
%   xy_c: center of the inclusions in each cell
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
%   M1: double (sparse), finite element stiffness matrix n1 * int_\partial D u*v dx
%
%   M2: double (sparse), finite element stiffness matrix n2 * int_\partial D u*v dx
% ------------------------------------------------------------------------------

n_nodes = size(p,1);
[n_elements, n_dof] = size(t);

% number of nodes on each inclusion boundary
N_inc0_left = size(xy_inc_left,1)/(n_sup*n_inc);
N_inc0_right = size(xy_inc_right,1)/(n_sup*n_inc);


% ==============================================================================
% boundary of inclusions
% nd: total number of nodes for each inclusion in each cell
% xy_inc: coordinates of each node for each inclusion in each cell
% map_inc: node number of each node for each inclusion in each cell
% ==============================================================================

nd = zeros(n_inc,2*n_sup);
xy_inc = zeros(1000,2*n_inc,2*n_sup);
map_inc = zeros(1000,n_inc,2*n_sup);
 
% label the inclusion boundary nodes
for n = 1:n_sup
    for l = 1:n_inc
        n_start_left = (n-1)*(N_inc0_left*n_inc) + (l-1)*N_inc0_left+1;
        n_end_left = (n-1)*(N_inc0_left*n_inc) + l*N_inc0_left;
        [~,label_on_left_inc] = inpolygon(p(:,1),p(:,2),xy_inc_left(n_start_left:n_end_left,1),xy_inc_left(n_start_left:n_end_left,2));
        idx_inc_bnd_nl = find(label_on_left_inc == 1);
        num_inc_nj = length(idx_inc_bnd_nl);
        
        hold on; plot(p(idx_inc_bnd_nl,1),p(idx_inc_bnd_nl,2),'ro');
        
        pi_loc = n;
        nd(l,pi_loc) = num_inc_nj;
        xy_inc(1:num_inc_nj,2*(l-1)+1:2*l,pi_loc) = p(idx_inc_bnd_nl,:);
        map_inc(1:num_inc_nj,l,pi_loc) = idx_inc_bnd_nl;
                 
    end
    for l = 1:n_inc
        n_start_right = (n-1)*(N_inc0_right*n_inc) + (l-1)*N_inc0_left+1;
        n_end_right = (n-1)*(N_inc0_right*n_inc) + l*N_inc0_left;
        [~,label_on_right_inc] = inpolygon(p(:,1),p(:,2),xy_inc_right(n_start_right:n_end_right,1),xy_inc_right(n_start_right:n_end_right,2));
         
        idx_inc_bnd_nl = find(label_on_right_inc == 1);
        num_inc_nj = length(idx_inc_bnd_nl);
         
        hold on; plot(p(idx_inc_bnd_nl,1),p(idx_inc_bnd_nl,2),'ro');

        pi_loc = n_sup + n;
        nd(l,pi_loc) = num_inc_nj;
        xy_inc(1:num_inc_nj,2*(l-1)+1:2*l,pi_loc) = p(idx_inc_bnd_nl,:);
        map_inc(1:num_inc_nj,l,pi_loc) = idx_inc_bnd_nl;
         
    end
end
 
 

% ==============================================================================
% Label nodes
% ==============================================================================
tol = 1e-4;                % tolerance should be smaller than mesh width
idx = zeros(n_nodes,1);
n_equations = 0;
 
x_min = xy_box(1); x_max = xy_box(2); 
y_min = xy_box(3); y_max = xy_box(4);
lat_c = y_max - y_min;

% label the nodes intially, not including nodes on the boundaries on the right
 
for i = 1:n_nodes     
    if ( ( abs(p(i,1)-x_min)<tol ) | ( abs(p(i,1)-x_max)<tol ) )
        % left or right bnd nodes: 0 (Dirichelt bnd cond)
        idx(i) = 0;
        % hold on; plot(p(i,1),p(i,2),'bo');
    elseif ( abs(p(i,2)-y_max)<tol )
        % top bnd nodes: -1
        idx(i) = -1;
        % hold on; plot(p(i,1),p(i,2),'ro');
    else
        % interior nodes and nodes on the inclusion bnd
        n_equations = n_equations + 1;
        idx(i) = n_equations;                                                  
    end
end

 
% sort the nodes on the inclusion boundary anti-clockwisely
for nc = 1:2*n_sup
    
    if (nc <= n_sup )
        shift_x = -nc*lat_c;
    else
        shift_x = (nc-n_sup-1)*lat_c;
    end
    
    for l = 1:n_inc
        xy_incl = [ xy_inc(1:nd(l,nc),2*l-1,nc) xy_inc(1:nd(l,nc),2*l,nc) ];
        zl = (xy_incl(:,1)-shift_x-xy_c(l,1)) + 1i*(xy_incl(:,2)-xy_c(l,2));
        thetal = angle(zl); 
        [~,Il] = sort(thetal);
        xy_inc(1:nd(l,nc),2*(l-1)+1:2*l,nc) = xy_incl(Il,:);
        map_inc(1:nd(l,nc),l,nc) = map_inc(Il,l,nc);
    end    
end


% periodic boundary mapping: map the nodes on top boundary to 
% the nodes on the bottom boundary

bot_node = find(abs(p(:,2)-y_min)<tol & abs(p(:,1)-x_min)>tol & abs(p(:,1)-x_max)>tol);
idx_bot = idx(bot_node);

for i = 1:n_nodes
    if idx(i) == -1      
        idx(i) = idx_bot( abs(p(i,1)-p(bot_node,1))<tol );
        % hold on; plot(p(i,1),p(i,2),'co');
        % temp = bot_node(abs(p(i,1)-p(bot_node,1))<tol);
        % hold on; plot(p(temp,1),p(temp,2),'co');
    end
end

if ( ~isempty(find(idx<0)) )
    disp('periodic boundary mapping not successful, adjust mesh size');
    return
end




% ==============================================================================
% Assembly in the elements
% ==============================================================================

n_nonzeros = n_dof^2*n_elements;  % estimate number of non-zero entries

mat  = zeros(n_nonzeros,1);    % ( dphi_i, dphi_j )
matx = zeros(n_nonzeros,1);    % ( phi_i, (mu*dphi_j)_x )
maty = zeros(n_nonzeros,1);    % ( phi_i, (mu*dphi_j)_x )
matu = zeros(n_nonzeros,1);    % ( phi_i, phi_j )
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



% ==============================================================================
% Assembly over the inclusions boundary
% ==============================================================================

mat_n1 = zeros(4*sum(nd),1);
mat_n2 = zeros(4*sum(nd),1);

ir_e = zeros(4*sum(nd),1);
ic_e = zeros(4*sum(nd),1);

count_e = 0;


% Assembly over the edges on inclusions
for nc = 1:2*n_sup
    for l = 1:n_inc
        for ne = 1:nd(l,nc)
    
            j1 = ne; j2 = (ne+1)*(ne<nd(l,nc)) + 1*(ne==nd(l,nc));
            xyj1 = [ xy_inc(j1,2*l-1,nc), xy_inc(j1,2*l,nc) ]; 
            xyj2 = [ xy_inc(j2,2*l-1,nc), xy_inc(j2,2*l,nc) ];
            hj = norm(xyj2-xyj1,2);
            vec = (xyj2-xyj1)/hj;
            nml = [ vec(2), -vec(1) ];
    
            hold on; plot( [xyj1(1) xyj2(1)], [xyj1(2) xyj2(2)], 'k-o' );
            quiver((xyj1(1)+xyj2(1))/2, (xyj1(2)+xyj2(2))/2, nml(1)/20, nml(2)/20, 'r' );
    
            mat_n1_loc = nml(1) * hj * [ 1/3 1/6; 1/6 1/3 ];
            mat_n2_loc = nml(2) * hj * [ 1/3 1/6; 1/6 1/3 ];
    
            for ni = 1:2
                for nj = 1:2            
                    count_e = count_e + 1;
                    ir_e(count_e) = idx(map_inc( j1*(ni==1)+j2*(ni==2), l, nc ));
                    ic_e(count_e) = idx(map_inc( j1*(nj==1)+j2*(nj==2), l, nc ));
            
                    mat_n1(count_e) = mat_n1_loc(ni,nj);
                    mat_n2(count_e) = mat_n2_loc(ni,nj); 
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

ir_e = ir_e(1:count_e);
ic_e = ic_e(1:count_e);
mat_n1 = mat_n1(1:count_e);
mat_n2 = mat_n2(1:count_e);
   
M1 = sparse(ir_e,ic_e,mat_n1,n_equations,n_equations);
M2 = sparse(ir_e,ic_e,mat_n2,n_equations,n_equations);
