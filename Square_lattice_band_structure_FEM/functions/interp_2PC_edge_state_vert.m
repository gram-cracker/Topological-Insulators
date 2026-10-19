%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                            %
%   Interpolate the Bloch mode at rectangular grid points    %      
%   for the joint structure                                  %                                                     
%                                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

addpath('./functions');

p=dlmread('../p.txt');
t=dlmread('../t.txt');
idx=dlmread('../idx.txt');
edge_state_n=dlmread('../edge_state_n.txt');

u = edge_state_n(:,1);
% twod_plotc(7,p,t,real(u));

% -------------------------------------------------------------------------
% Geometry description:
% -------------------------------------------------------------------------

% supercell vertices

n_sup = 8;                           % number of cells in top/bottom structures
lat_c = 1;                           % square lattice constant 
Lx = lat_c;  Ly = n_sup * lat_c;     % supercell size for the top/bottom structures

xy1 =  [0,-Ly];
xy2 =  [2*Lx,-Ly];
xy3 =  [2*Lx,Ly];
xy4 =  [0,Ly];

cell_v = [ xy1; xy2; xy3; xy4 ];

x_corner_left = cell_v(1,:);
x_corner_right = cell_v(3,:);

% rectangular grid points over the domain
Mx = 50; My = 400;
x_vec = linspace(x_corner_left(1),x_corner_right(1),Mx);
y_vec = linspace(x_corner_left(2),x_corner_right(2),My);

[XX,YY] = meshgrid(x_vec,y_vec);
bloch_mode = zeros(size(XX));


% Bloch wave vector
kappa_parallel = 0.0 * 2*pi/lat_c * [1, 0];

% -------------------------------------------------------------------------
% interpolate the solution to obtain the Bloch mode
% -------------------------------------------------------------------------

for m = 1:Mx
    for n = 1:My
              
        m        
        x = XX(n,m); y = YY(n,m);
                    
        % map (x,y) to (x_bar,y_bar) in the reference period
        x_bar = mod(x,Lx);
        y_bar = y;
                                
        bloch_mode(n,m) = exp(1i*kappa_parallel*[x,y]') * u_interp(u, p, t, x_bar, y_bar);
    
    
    end 
end
% -------------------------------------------------------------------------

figure;  pcolor(XX,YY,real(bloch_mode))
shading interp; axis equal; axis tight; axis off; 
caxis([-1, 1]); colormap jet;