% -------------------------------------
% square inclusions in the cells
% -------------------------------------

n_inc = 1;                         % number of inclusions
L_inc = 0.45;                      % size of the inclusion

theta = pi/18;                     % rotation angle

% rotation matrix on the left

% inclusion in the first left cell
xy_inc0_left = [ lat_c/2, lat_c/2 ]' + [ -L_inc/2, L_inc/2, -L_inc/2; -L_inc/2, -L_inc/2, L_inc/2 ];
xy_inc0_left = xy_inc0_left';

% all the inclusions in the left super cells
xy_inc_left = zeros(3*n_sup,2);
for n = 1:n_sup
    xy_inc_left(3*(n-1)+1:3*n,:) = [-n*lat_c, 0] + xy_inc0_left; 
end

% inclusion in the first right cell
xy_inc0_right = [ lat_c/2, lat_c/2 ]' + [ L_inc/2, L_inc/2, -L_inc/2; -L_inc/2, L_inc/2, L_inc/2 ];
xy_inc0_right = xy_inc0_right';

% all the inclusions in the left super cells
xy_inc_right = zeros(3*n_sup,2);
for n = 1:n_sup
    xy_inc_right(3*(n-1)+1:3*n,:) = [ (n-1)*lat_c, 0] + xy_inc0_right;
end
