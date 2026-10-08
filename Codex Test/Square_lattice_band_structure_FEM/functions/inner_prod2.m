function val = inner_prod2(p,t,t_num,eps_inclu,eps0,pol,u1,u2)
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
    u1_loc = u1(nodes_loc);
    u2_loc = u2(nodes_loc);
    
    % permittivity
    if ( pol == 'TE' )
        flag_inc = (t_num(ne) == 1);
        eps_loc = flag_inc * eps_inclu + (~flag_inc) * eps0;
    elseif ( pol == 'TM' )
        eps_loc = 1;
    else
        display('wrong polarization');
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
    
    val = val + eps_loc * u2_loc.' * (matu_loc2 * u1_loc);

end

end

