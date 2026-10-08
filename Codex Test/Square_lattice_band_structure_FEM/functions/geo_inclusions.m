

%%%%%%%%%%%%%%%%%%%%%%%%%
% circular inclusions 1 %
%%%%%%%%%%%%%%%%%%%%%%%%%

% circular inclusion
% n_inc = 1;                   % number of inclusions
% R = 0.2;                    % radius of the inclusion
% xy_c = 1/2 * [ Lx; Ly ];
% 
% ntheta = 50;
% dtheta = 2*pi/ntheta;
% theta = 0:dtheta:2*pi-dtheta;
% 
% xy_inc = zeros(n_inc*ntheta,2);
% for j = 1:n_inc
%     xy_inc((j-1)*ntheta+1:j*ntheta,:) = [  xy_c(1) + R * cos(theta); xy_c(2) + R * sin(theta) ]';
% end



%%%%%%%%%%%%%%%%%%%%%%%%%
% circular inclusions 2 %
%%%%%%%%%%%%%%%%%%%%%%%%%

% rot_ang = -pi/18;
% 
% % inclusions in the corners
% n_inc1 = 4;                              % number of inclusions in the corners
% R1 = 0.05 * lat_c;                       % radius of the inclusions
% d1 = sqrt(2) * lat_c/4;                  % distance between the orgin and the inclusion center
% theta_c1 = rot_ang + pi/4 + 2*pi/n_inc1 * ( 0:1:n_inc1-1 );
% xy_c1 = lat_c/2 + d1 * [ cos(theta_c1); sin(theta_c1) ]';
% 
% % inclusions on the sides
% n_inc2 = 4;                              % number of inclusions on the sides
% R2 = 0.075 * lat_c;                      % radius of the inclusions
% d2 = lat_c/4;                            % distance between the orgin and the inclusion center
% theta_c2 = rot_ang + 2*pi/n_inc2 * ( 0:1:n_inc2-1 );
% xy_c2 = lat_c/2 + d2 * [ cos(theta_c2); sin(theta_c2) ]';
% 
% ntheta = 50;
% dtheta = 2*pi/ntheta;
% theta = 0:dtheta:2*pi-dtheta;
% xy_inc1 = zeros(n_inc1*ntheta,2);
% xy_inc2 = zeros(n_inc2*ntheta,2);
% 
% 
% for j = 1:n_inc1
%     xy_inc1((j-1)*ntheta+1:j*ntheta,:) = [ xy_c1(j,1) + R1 * cos(theta);  xy_c1(j,2) + R1 * sin(theta) ]';
%     xy_inc2((j-1)*ntheta+1:j*ntheta,:) = [ xy_c2(j,1) + R2 * cos(theta);  xy_c2(j,2) + R2 * sin(theta) ]';
%     % hold on; plot(xy_inc1((j-1)*ntheta+1:j*ntheta,1), xy_inc1((j-1)*ntheta+1:j*ntheta,2) );
%     % hold on; plot(xy_inc2((j-1)*ntheta+1:j*ntheta,1), xy_inc2((j-1)*ntheta+1:j*ntheta,2) );
% end
% 
% n_inc = n_inc1 + n_inc2;
% xy_inc = [ xy_inc1; xy_inc2 ];


%%%%%%%%%%%%%%%%%%%%%%%%%
% circular inclusions 3 %
%%%%%%%%%%%%%%%%%%%%%%%%%

% rot_ang = 0; % -pi/8;
% 
% % inclusions in the corners
% n_inc = 4;                              % number of inclusions in the corners
% R = 0.12 * lat_c;                       % radius of the inclusions
% d = sqrt(2) * lat_c/4;                  % distance between the orgin and the inclusion center
% theta_c = rot_ang + pi/4 + 2*pi/n_inc * ( 0:1:n_inc-1 );
% xy_c = lat_c/2 + d * [ cos(theta_c); sin(theta_c) ]';
% 
% ntheta = 50;
% dtheta = 2*pi/ntheta;
% theta = 0:dtheta:2*pi-dtheta;
% xy_inc = zeros(n_inc*ntheta,2);
% 
% 
% for j = 1:n_inc
%     xy_inc((j-1)*ntheta+1:j*ntheta,:) = [ xy_c(j,1) + R * cos(theta);  xy_c(j,2) + R * sin(theta) ]';
% end



%%%%%%%%%%%%%%%%%%%%%%%%%
% circular inclusions 4 %
%%%%%%%%%%%%%%%%%%%%%%%%%

% inclusions in the corners
n_inc = 2;                   % number of inclusions
R = 0.12 * lat_c;            % radius of the inclusions
d = 0.12 * lat_c; % 0.12; 0.35 % distance between the cell center and the inclusion center
xy_c = zeros(2,2);
xy_c(1,:) = [ Lx/2; Ly/2 ] + [ d; d ];
xy_c(2,:) = [ Lx/2; Ly/2 ] - [ d; d ];

ntheta = 50;
dtheta = 2*pi/ntheta;
theta = 0:dtheta:2*pi-dtheta;
xy_inc = zeros(n_inc*ntheta,2);

for j = 1:n_inc
    xy_inc((j-1)*ntheta+1:j*ntheta,:) = [ xy_c(j,1) + R * cos(theta'),  xy_c(j,2) + R * sin(theta') ];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% square/rectangular inclusions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n_inc = 1;
% Lx_inc = 0.45; Ly_inc = 0.25; % 0.25;
% 
% dtheta = 0; % -(pi/2+pi/8); % pi/8;
% theta = pi/4 + dtheta;
% Rot = [  cos(theta) sin(theta);
%         -sin(theta) cos(theta) ];
% 
% shift = Ly_inc * sqrt(2)/8 * [ 1, 1 ]';    
% xy_inc = [ Lx/2, Ly/2 ]' + shift + Rot * [ -Lx_inc/2, Lx_inc/2, Lx_inc/2, -Lx_inc/2; -Ly_inc/2, -Ly_inc/2, Ly_inc/2 Ly_inc/2 ];
% xy_inc = xy_inc';


%%%%%%%%%%%%%%%%%%%%%%%%%
% triangular inclusions %
%%%%%%%%%%%%%%%%%%%%%%%%%
% n_inc = 1;
% Lx_inc = 0.45; Ly_inc = 0.45;
% 
% % triangle pointing up
% xy_inc = [ Lx/2, Ly/2 ]' + [ -Lx_inc/2, Lx_inc/2, -Lx_inc/2; -Ly_inc/2, -Ly_inc/2, Ly_inc/2 ];
% 
% % tiranlge pointing down
% % xy_inc = [ Lx/2, Ly/2 ]' + [ Lx_inc/2, Lx_inc/2, -Lx_inc/2; -Ly_inc/2, Ly_inc/2, Ly_inc/2 ];
% 
% xy_inc = xy_inc';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% example: optimal structures with circular inclusion as initial guess %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load('./data/opt_joint_PC_opt_circ_inc.mat');
% V = reshape(V2_saved{10},64,64);       % data for two PCs: V1 and V2
% Vinterp = scatteredInterpolant(X(:),Y(:),V(:),'nearest');
% Vinterp2 = @(x,y) Vinterp( mod(Minv(1,1)*x + Minv(1,2)*y,1), mod(Minv(2,1)*x + Minv(2,2)*y,1) );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% example: optimal structures with triangular inclusion as initial guess %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load('./data/opt_joint_PC_opt_tri_inc_Berry_constraint_single_k.mat');
% % load('./data/opt_joint_PC_init_guess_tri_inc.mat');
% V = reshape(V2,64,64);     % data for two PCs: V1 and V2
% Vinterp = scatteredInterpolant(X(:),Y(:),V(:),'nearest');
% Vinterp2 = @(x,y) Vinterp( mod(Minv(1,1)*x + Minv(1,2)*y,1), mod(Minv(2,1)*x + Minv(2,2)*y,1) );


