function val = inner_prod_interp(p,t,t_num,pol,u1,u2,Vinterp2)
%INNER_PROD computes the inner product of two functions u1 and u2


n_nodes = size(p,1);

if (length(u1)~=length(u2) | length(u1) ~= n_nodes)
    display('the length of u1 and u2 is wrong');
    return;
end
    
[n_elements, n_dof] = size(t);

val = 0;

[r,s,w] = twod_gauss(7);       % Gauss nodes and weights on standard element

for ne = 1:n_elements
    
    nodes_loc = t(ne,:);
    xy_loc = p(nodes_loc,:);
    xy_center = sum(p(nodes_loc,:))/3;
    u1_loc = u1(nodes_loc);
    u2_loc = u2(nodes_loc);
    
    % permittivity
    if ( pol == 'TE')       
        eps_loc = Vinterp2(xy_center(1),xy_center(2));
    elseif ( pol == 'TM')
        eps_loc = 1;
    else
        disp('Wrong EM ploarization');
        return;
    end
    
    % Compute local gauss points/weights and shape functions there
    [~,wg,phi,~,~] = twod_shape(xy_loc,r,s,w);
    
    matu_loc = twod_bilinear(eps_loc,phi,phi,wg);
    
    xi = p(nodes_loc(1),1); yi = p(nodes_loc(1),2); 
    xj = p(nodes_loc(2),1); yj = p(nodes_loc(2),2); 
    xk = p(nodes_loc(3),1); yk = p(nodes_loc(3),2); 
 
    area_e=det([1 xi yi; 1 xj yj; 1 xk yk])/2;
    matu_loc2 = area_e/12 * [ 2 1 1; 1 2 1; 1 1 2 ];
    
    val = val + eps_loc * u2_loc' * (matu_loc2 * u1_loc);
    % val = val + eps_loc * u1_loc' * (matu_loc2 * u2_loc);
    
end

end

