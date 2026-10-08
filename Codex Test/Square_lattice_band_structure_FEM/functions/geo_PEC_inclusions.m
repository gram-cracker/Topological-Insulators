
% square inclusion
% lx_inc = 0.35; ly_inc = 0.35;
% % lx_inc = 0.35; ly_inc = 0.35;
% hhx = 0.01; hhy = 0.01;     % (smaller than the mesh size hmax)
% NNx = lx_inc/hhx; NNy = ly_inc/hhy;
% xy_inc = [  zeros(1,NNy),     hhx:hhx:lx_inc,  lx_inc*ones(1,NNy), lx_inc-hhx:-hhx:0; ...
%             ly_inc-hhy:-hhy:0, zeros(1,NNx),    hhy:hhy:ly_inc,        ly_inc*ones(1,NNx), ]';
% xy_inc(:,1) = (Lx/2-lx_inc/2) + xy_inc(:,1);        
% xy_inc(:,2) = (Ly/2-ly_inc/2) + xy_inc(:,2);
% 
% xy_center = [ Lx/2; Ly/2 ];
% 
% % rotation of the inclusion
% theta = 0; % pi/6;
% R = [ cos(theta) sin(theta);
%      - sin(theta)  cos(theta) ];
% xy_inc = ( R * (xy_inc' - xy_center) + xy_center )';


% circular inclusion
xc = [ Lx/2 Ly/2 ]; 
R = 0.1;
dtheta = 2*pi/100;
theta = 0:dtheta:2*pi-dtheta;
xy_inc = [ xc(1) + R * cos(theta);  xc(2) + R * sin(theta) ]';

% plot(xy_cell(:,1),xy_cell(:,2));
% hold on;
% plot(xy_inc(:,1),xy_inc(:,2));
