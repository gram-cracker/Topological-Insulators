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
Lx = n_sup * lat_c; Ly = 2*lat_c;    % supercell size for the top/bottom structures

xy1 =  [-Lx,0];
xy2 =  [Lx,0];
xy3 =  [Lx,Ly];
xy4 =  [0,Ly];

cell_v = [ xy1; xy2; xy3; xy4 ];

x_corner_left = cell_v(1,:);
x_corner_right = cell_v(3,:);

% rectangular grid points over the domain
Mx = 400; My = 50;
x_vec = linspace(x_corner_left(1),x_corner_right(1),Mx);
y_vec = linspace(x_corner_left(2),x_corner_right(2),My);

[XX,YY] = meshgrid(x_vec,y_vec);
bloch_mode = zeros(size(XX));


% Bloch wave vector
kappa_parallel = 0.0 * 2*pi/lat_c * [0, 1];

% -------------------------------------------------------------------------
% interpolate the solution to obtain the Bloch mode
% -------------------------------------------------------------------------

warning('off','all');
parpool('local',28);

parfor m = 1:Mx
    for n = 1:My
              
        m        
        x = XX(n,m); y = YY(n,m);
                    
        % map (x,y) to (x_bar,y_bar) in the reference period
        x_bar = x;
        y_bar = mod(y,lat_c);
                                
        bloch_mode(n,m) = exp(1i*kappa_parallel*[x,y]') * u_interp(u, p, t, x_bar, y_bar);
    
    
    end 
end
% -------------------------------------------------------------------------

poolobj = gcp('nocreate');
delete(poolobj);

figure;  pcolor(XX,YY,real(bloch_mode))
title('$k_\parallel=0.0$, $\omega=0.5962$')
shading interp; axis equal; axis tight; %axis off; 
caxis([-1, 1]); colormap jet;

dlmwrite('../XX.txt',XX); dlmwrite('../YY.txt',YY);
dlmwrite('../Bloch_mode2_kappa_0.3.txt',bloch_mode);