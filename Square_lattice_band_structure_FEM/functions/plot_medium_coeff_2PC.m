

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                %
%   plot medium coefficient for the joint PC     %      
%                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

% -------------------------------------------------------------------------
% Geometry description:
% -------------------------------------------------------------------------

lat_c = 1;  % lattice constant
n_sup = 8;  % number of super cells
n_py = 2;   % number of periods along the y direction
lx = lat_c*n_sup;
ly = lat_c*n_py;

% supercell vertices
xy1 =  [-lx,0];
xy2 =  [lx,0];
xy3 =  [lx,ly];
xy4 =  [-lx,ly];

cell_v = [ xy1; xy2; xy3; xy4 ];

x_corner_left = cell_v(1,:);
x_corner_right = cell_v(3,:);

Mx = 2000; My = 200;
x_vec = linspace(x_corner_left(1),x_corner_right(1),Mx);
y_vec = linspace(x_corner_left(2),x_corner_right(2),My);

[XX,YY] = meshgrid(x_vec,y_vec);
C_matrix = zeros(size(XX));

% box for the joint structure
% xy10 = [0 0] - [lat_c * n_sup, 0];
% xy20 = [0 0] + [lat_c * n_sup, 0];
% xy30 = [lat_c * n_sup, 0] * [ 1/2, sqrt(3)/2 ];
% xy40 = n_py * lat_c * [1/2 sqrt(3)/2]
% 
% parallel_box = [ xy10; xy20; xy30; xy40 ];
% figure; fill(parallel_box(:,1),parallel_box(:,2),'c','Linestyle','no');

% slope = sqrt(3);
% a0 = sqrt(2/sqrt(3)); % lattice constant in the imported example
% scal = a0/lat_c;

% Example 4
% -------------------------------------------------------------------------
load('../data/opt_joint_PC_opt_tri_inc_Berry_constraint_single_k.mat');
Vinterp1 = scatteredInterpolant(X(:),Y(:),V1(:),'nearest');
Vinterp2 = scatteredInterpolant(X(:),Y(:),V2(:),'nearest');
shift1 =  [ 0, 0 ]; 
shift2 =  [ 0, 0 ];
% -------------------------------------------------------------------------

for m = 1:Mx
    for n = 1:My
               
        x = XX(n,m); y = YY(n,m);
                
        flag1 = x>0;
    
        if ( flag1 == 1 )
            Vinterp = Vinterp2;
            shift = shift1;
        else
            Vinterp = Vinterp1;
            shift = shift2;
        end
        
        Vinterp_fun = @(x,y) Vinterp( mod(Minv(1,1)*x + Minv(1,2)*y,1), mod(Minv(2,1)*x + Minv(2,2)*y,1) );
        coeff_c_loc = Vinterp_fun(x+shift(1),y+shift(2));
        
        C_matrix(n,m) = coeff_c_loc;
            
        % if (abs(coeff_c_loc- 11.7) < 1e-3 )
        % hold on; plot(x,y,'.');
        % end                  
    
    end 
end

% axis equal
% axis tight
% axis off

figure;  pcolor(XX,YY,C_matrix)
shading interp; axis equal; axis tight; axis off; % colormap bone;