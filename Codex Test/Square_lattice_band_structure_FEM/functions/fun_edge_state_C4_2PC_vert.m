function [omega_edge, U_edge] = fun_edge_state_C4_2PC_vert(omega,gap_min,gap_max,U,p,t,idx,Ly)
% fun_edge_state identifies the edge mode among all the computed
% eigenfunctions

n_nodes = size(p,1);
N = size(U,2);
eps = 0.35;

j = 0;
omega_edge = zeros(size(omega));
U_edge = zeros(n_nodes,N);

for n = 1:N
    
    U(:,n) = U(:,n)/max(abs(U(:,n)));    

    un = zeros(n_nodes,1);
    for i = 1:n_nodes
        if (idx(i)>0)
            un(i) = U(idx(i),n);
        end
    end

    un_norm = norm(un,2);
    
    idx_nbd = find(abs(p(:,2))>Ly/2);
    un2 = un(idx_nbd);
    un_norm2 = norm(un2,2); 
    
    if ((un_norm2/un_norm<eps) & (omega(n)>gap_min) & (omega(n)<gap_max))
        
        j = j+1;
        omega_edge(j) = omega(n);
        U_edge(:,j) = un;
                
        % twod_plotc(j,p,t,(-1)^(j+1)*real(un));
        % title([ '$k_\parallel = 0.5$, ', '$\omega$ = ', num2str(real(omega(n))) ]);
        % axis tight;
        
    end
    
end

omega_edge = omega_edge(1:j);
U_edge = U_edge(:,1:j);

end

