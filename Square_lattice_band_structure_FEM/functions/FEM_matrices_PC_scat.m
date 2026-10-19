function [M,Me,An,F,idx,l2r] = FEM_matrices_PC_scat(p,t,t_num,neu_nodes,xy_box,eps_inc,eps_bk,mu_inc,mu_bk,pol,pml_para,k0,theta)
% Description: compute the mass and stiffness matrices
%
% Inputs:
%
%   p: double, (2,n_nodes) array of mesh nodes
%
%   t: int, (n_elements, n_dof) triangular element matrix
%
%   t_num: array of regions numbers for each triangle in t
%
%   neu_nodes: Neumann nodes 
%
%   xy_box: double, [x_min, x_max, y_min, y_max] denote the domain
%
%   eps_inc/mu_inc: permittivity and permeability of the inclusions
%
%   eps_bk/mu_bk: permittivity and permeability of the background
%
%   pol: TE/TM polarization
%
%   pml_para = [ y0, ly]: starting point and thickness of pml
%
%   k0: free-space wavenumber
%
%   theta: incident angle
%
% Outputs:
%
%   M: double, (sparse), finite element mass matrix int_D a * u'*v' dx
%
%   Me: double (sparse), finite element stiffness matrix int_D c * u'*v' dx
%
%   An: double,  int_Ie phi_j(x) exp(i*kx*x) dx
%
%   Fn: double, int_[0,L] phi_j(x) exp(-i*kx(n)*x) dx
%
% ------------------------------------------------------------------------------

k_top = k0*sqrt(eps_bk(1)*mu_bk(1));
k_bot = k0*sqrt(eps_bk(2)*mu_bk(2));

kx = k0*sin(theta);
n_nodes = size(p,1);
[n_elements, n_dof] = size(t);

x_min = xy_box(1); x_max = xy_box(2);
y_min = xy_box(3); y_max = xy_box(4);

Lx = x_max - x_min;

if (pol == 'TE')
    coeff_a = 1/mu_inc;   coeff_c = eps_inc;
    coeff_a0_top = 1/mu_bk(1);   coeff_c0_top = eps_bk(1);
    coeff_a0_bot = 1/mu_bk(2);   coeff_c0_bot = eps_bk(2);
elseif (pol == 'TM')
    coeff_a = 1/eps_inc;  coeff_c_top = mu_inc;
    coeff_a0_top = 1/eps_bk(1);  coeff_c0_top = mu_bk(1);
    coeff_a0_bot = 1/eps_bk(2);  coeff_c0_bot = mu_bk(2);
else
    display('Wrong EM ploarization');
    return;
end

% ==============================================================================
% Label nodes
% ==============================================================================
 tol = 1e-3;                % tolerance should be smaller than mesh width
 idx = zeros(n_nodes,1);
 n_equations = 0;
 
 % label the nodes intially, not including nodes on the right and top boundaries
 
 for i = 1:n_nodes
     if ( abs(p(i,2)-y_min)<tol )
         % bottom bnd nodes: 0 (Dirichelt bnd cond)
         idx(i) = 0;
     elseif ( (abs(p(i,1)-x_max)<tol) && (abs(p(i,2)-y_max)>tol) && (abs(p(i,2)-y_min)>tol) )
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
l2r = zeros(n_nodes,1);

for i = 1:n_nodes
    if idx(i) == -1
        l2r(i) = idx_left( abs(p(i,2)-p(left_node,2))<tol );
        % hold on; plot(p(i,1),p(i,2),'go');
        % hold on; plot(p(left_node,1),p(left_node,2),'go');
    end
end




% Assembly
% ==============================================================================

n_nonzeros = n_dof^2*n_elements;  % estimate number of non-zero entries

mat  = zeros(n_nonzeros,1);       % ( a*dphi_i, dphi_j )
mate = zeros(n_nonzeros,1);       % ( c*phi_i, phi_j )
ir = zeros(n_nonzeros,1);         % row indices
ic = zeros(n_nonzeros,1);         % column indices

[r,s,w] = twod_gauss(7);          % Gauss nodes and weights on standard element
count = 0;


for ne = 1:n_elements
    
    nodes_loc = t(ne,:);
    xy_loc = p(nodes_loc,:);
    xy_center = sum(xy_loc,1)/3;
    
    % Compute local gauss points/weights and shape functions there
    [~,wg,phi,px,py] = twod_shape(xy_loc,r,s,w);
        
    % permittivity and permeability
    flag_inc = (t_num(ne) == 1);
    flag_bk_top = (t_num(ne) == 2) & (xy_center(2)>0);
    flag_bk_bot = (t_num(ne) == 2) & (xy_center(2)<=0);

    coeff_a_loc = flag_inc * coeff_a + flag_bk_top * coeff_a0_top + flag_bk_bot * coeff_a0_bot;
    coeff_c_loc = flag_inc * coeff_c + flag_bk_top * coeff_c0_top + flag_bk_bot * coeff_c0_bot;
       
    if ( abs(coeff_c_loc - 1) < 1e-3 )
         hold on; plot(xy_center(1),xy_center(2),'ro');
    elseif ( abs(coeff_c_loc - 4) < 1e-3 )
         hold on; plot(xy_center(1),xy_center(2),'go');
    end

    % PML absorbing parameters
    sy = fun_sy(xy_center(2),pml_para(1),pml_para(2),k_bot);
    
    % local stiff matrix
    mat_loc = coeff_a_loc*sy*twod_bilinear(1,px,px,wg) + coeff_a_loc*(1/sy)*twod_bilinear(1,py,py,wg);
    mate_loc = coeff_c_loc*sy*twod_bilinear(1,phi,phi,wg);

    for ni = 1:n_dof
 
        idx_ni = idx(nodes_loc(ni));
        coeff_ni = 1.0;
        
        if (idx_ni == -1)
            % map to right boundary node to the left
            idx_ni = l2r(nodes_loc(ni));
            coeff_ni = exp(1i*kx*Lx);
        end
        
        if (idx_ni > 0)
            for nj = 1:n_dof

                idx_nj = idx(nodes_loc(nj));
                coeff_nj = 1.0;

                if (idx_nj == -1)
                    % map to right boundary node to the left
                    idx_nj = l2r(nodes_loc(nj));
                    coeff_nj = exp(1i*kx*Lx);
                end
                if (idx_nj > 0)
                    count = count + 1;
                    ir(count) = idx_ni;
                    ic(count) = idx_nj;
                    mat(count)  = conj(coeff_ni)*coeff_nj*mat_loc(ni,nj);
                    mate(count) = conj(coeff_ni)*coeff_nj*mate_loc(ni,nj);
                end
            end  % for nj
        end  % if idx_ni
    end  % for ni

    
end


% keep only nonzero entries
ir = ir(1:count);
ic = ic(1:count);
mat = mat(1:count);
mate = mate(1:count);

M = sparse(ir,ic,mat,n_equations,n_equations);
Me = sparse(ir,ic,mate,n_equations,n_equations);



% ------------------------------------------------------------------------------
% F(n,j) = int_[0,L] phi_j(x) exp(-i*k1(n)*x) dx for j on top edge, |n| <= n_max 
% ------------------------------------------------------------------------------

n_max = 2*floor(k0*Lx/(2*pi))+1;  % max number of modes kept
n  = (-n_max:n_max)';
kxn = kx+2*pi*n/Lx;
flag_n = find(abs(kxn)<eps);   idx_n = find(abs(kxn)>eps);

% Edge connectivity (specific to pw linear)
n_neu_nodes = length(neu_nodes);
neu_edge_counter = [1:n_neu_nodes-1;2:n_neu_nodes]';
[n_neu_edges,n_neu_dof] = size(neu_edge_counter);

F = zeros(2*n_max+1,n_neu_nodes);
An = zeros(n_neu_nodes,1);

for n_edge = 1:n_neu_edges
    
    % Get local nodes and shape functions
    nodes_count_loc = neu_edge_counter(n_edge,:);
    nodes_loc = neu_nodes(nodes_count_loc);
    x_loc = p(nodes_loc,1);
    x1 = x_loc(1); x2 = x_loc(2);
    l = abs(x2-x1);
    
    % Compute local An = Ie exp(ik_x*x)*phi(x) dx
    An_loc = zeros(1,2);
    if (kx~=0)
        s = kx*l;
        An_loc(2) = l*exp(1i*kx*x1).*((1/s^2-1i/s)*exp(1i*s)-1/s^2);
        An_loc(1) = l*exp(1i*kx*x1).*(1i/s*(1-exp(1i*s)))- An_loc(2);
    else
        An_loc = l/2*[1 1];
    end
        
    % Local F = Ie phi(x) exp(1i*kxn*x) dx
    F_loc = zeros(2*n_max+1,2);
    
    s = -kxn*l;
    if (~isempty(flag_n))
        F_loc(flag_n,:) = l/2*[1 1];
    end
    if (~isempty(idx_n))
        F_loc(idx_n,2) = l*exp(-1i*kxn(idx_n)*x1).*((1./s(idx_n).^2-1i./s(idx_n)).*exp(1i*s(idx_n))-1./s(idx_n).^2);
        F_loc(idx_n,1) = l*exp(-1i*kxn(idx_n)*x1).*(1i./s(idx_n).*(1-exp(1i*s(idx_n))))- F_loc(idx_n,2);
    end
  
    % Local to global mapping
    for nt = 1:n_neu_dof
        n_col = nodes_count_loc(nt);
        F(:,n_col) = F(:,n_col) + F_loc(:,nt);
        An(n_col) = An(n_col) + An_loc(nt);
    end

end