function F = FEM_src_vec_C4(p,t,xy_box,xy_src,sign_src)
% Description: compute the mass and stiffness matrics
%
% Inputs:
%
%   p: double, (2,n_nodes) array of mesh nodes
%
%   t: int, (n_elements, n_dof) triangular element matrix
%
%   xy_box: double, [x_min, x_max, y_min, y_max] denote the domain
%
%   xy_src, sign_src: source locations and strength
%
% Outputs:
%
%   F: double, finite element source vector int_D f*v dx
% ------------------------------------------------------------------------------

% global xy_src sign_src

num_src = length(sign_src);
n_nodes = size(p,1);

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

F = zeros(n_equations,1);

for n = 1:num_src
    
    % Source function
    ne_src = fun_xy_loc( xy_src(n,1), xy_src(n,2), p, t );
    nodes_loc = t(ne_src,:);
    xi = p(nodes_loc(1),1); yi = p(nodes_loc(1),2); 
    xj = p(nodes_loc(2),1); yj = p(nodes_loc(2),2); 
    xk = p(nodes_loc(3),1); yk = p(nodes_loc(3),2); 
 
    area_e=det([1 xi yi; 1 xj yj; 1 xk yk])/2;

    xy_center = 1/3 * [ xi+xj+xk, yi+yj+yk ];
    % plot(xy_center(1),xy_center(2),'ro');
    
    % u=ui*Ni(x,y)+uj*Nj(x,y)+uk*Nk(x,y), where Ni(x,y)=(ai+bi*x+ci*y)/2area_e; 
    ai=xj*yk-xk*yj; bi=yj-yk; ci=xk-xj;
    aj=xk*yi-xi*yk; bj=yk-yi; cj=xi-xk;
    ak=xi*yj-xj*yi; bk=yi-yj; ck=xj-xi;
           
    Ni = (ai+bi*xy_src(n,1)+ci*xy_src(n,2))/(2*area_e);
    Nj = (aj+bj*xy_src(n,1)+cj*xy_src(n,2))/(2*area_e);
    Nk = (ak+bk*xy_src(n,1)+ck*xy_src(n,2))/(2*area_e);     

    F(idx(nodes_loc(1:3))) = F(idx(nodes_loc(1:3))) + sign_src(n) * [ Ni, Nj, Nk ]';


end
