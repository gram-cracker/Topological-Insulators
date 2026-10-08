% -------------------------------------
% square inclusions in the cells
% -------------------------------------

n_inc = 1;                         % number of inclusions
L_inc = 0.35;                       % size of the inclusion

xy_center = [ lat_c/2, lat_c/2 ];

theta = pi/6;                     % rotation angle

% rotation matrix on the left
Rot_left = [  cos(-theta) sin(-theta);    
             -sin(-theta) cos(-theta) ];

% inclusion in the first left cell
xy_inc0_left = xy_center' + Rot_left * [ -L_inc/2, L_inc/2, L_inc/2, -L_inc/2; -L_inc/2, -L_inc/2, L_inc/2 L_inc/2 ];
xy_inc0_left = xy_inc0_left';

% all the inclusions in the left super cells
xy_inc_left = zeros(4*n_sup,2);
for n = 1:n_sup
    xy_inc_left(4*(n-1)+1:4*n,:) = [-n*lat_c, 0] + xy_inc0_left; 
end


Rot_right = [  cos(theta) sin(theta);
              -sin(theta) cos(theta) ];

% inclusion in the first right cell
xy_inc0_right = xy_center' + Rot_right * [ -L_inc/2, L_inc/2, L_inc/2, -L_inc/2; -L_inc/2, -L_inc/2, L_inc/2 L_inc/2 ];
xy_inc0_right = xy_inc0_right';

% all the inclusions in the left super cells
xy_inc_right = zeros(4*n_sup,2);
for n = 1:n_sup
    xy_inc_right(4*(n-1)+1:4*n,:) = [ (n-1)*lat_c, 0] + xy_inc0_right;
end
