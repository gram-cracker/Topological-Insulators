

% ---------------------------------
% circular inclusions in the cells
% ---------------------------------


% inclusions in the corners
n_inc1 = 4;                              % number of inclusions in the corners
R1 = 0.05 * lat_c;                       % radius of the inclusions
d1 = sqrt(2) * lat_c/4;                  % distance between the orgin and the inclusion center

% inclusions on the sides
n_inc2 = 4;                              % number of inclusions on the sides
R2 = 0.075 * lat_c;                      % radius of the inclusions
d2 = lat_c/4;                            % distance between the orgin and the inclusion center

n_inc = n_inc1 + n_inc2;
ntheta = 30;
dtheta = 2*pi/ntheta;
theta = 0:dtheta:2*pi-dtheta;


% ---------------------------------------------
% inclusions on the left super cells
% ---------------------------------------------

% inclusion in the first left cell
rot_ang_left = pi/18;      % rotation angle
theta_c1_left = rot_ang_left + pi/4 + 2*pi/n_inc1 * ( 0:1:n_inc1-1 );
xy_c1_left = lat_c/2 + d1 * [ cos(theta_c1_left); sin(theta_c1_left) ]';

theta_c2_left = rot_ang_left + 2*pi/n_inc2 * ( 0:1:n_inc2-1 );
xy_c2_left = lat_c/2 + d2 * [ cos(theta_c2_left); sin(theta_c2_left) ]';

xy_inc1_left0 = zeros(n_inc1*ntheta,2);
xy_inc2_left0 = zeros(n_inc2*ntheta,2);
for j = 1:n_inc1
    xy_inc1_left0((j-1)*ntheta+1:j*ntheta,:) = [ xy_c1_left(j,1) + R1 * cos(theta);  xy_c1_left(j,2) + R1 * sin(theta) ]';
    xy_inc2_left0((j-1)*ntheta+1:j*ntheta,:) = [ xy_c2_left(j,1) + R2 * cos(theta);  xy_c2_left(j,2) + R2 * sin(theta) ]';
end

xy_inc_left0 = [ xy_inc1_left0; xy_inc2_left0 ];
len = size(xy_inc_left0,1);

% all the inclusions in the left super cells
xy_inc_left = zeros(n_sup*len,2);
for n = 1:n_sup
    xy_inc_left((n-1)*len+1:n*len,:) = [-n*lat_c, 0] + xy_inc_left0;
end


% ---------------------------------------------
% inclusions on the right super cells
% ---------------------------------------------

% inclusion in the first right cell
rot_ang_right = -pi/18;          % rotation angle
theta_c1_right = rot_ang_right + pi/4 + 2*pi/n_inc1 * ( 0:1:n_inc1-1 );
xy_c1_right = lat_c/2 + d1 * [ cos(theta_c1_right); sin(theta_c1_right) ]';

theta_c2_right = rot_ang_right + 2*pi/n_inc2 * ( 0:1:n_inc2-1 );
xy_c2_right = lat_c/2 + d2 * [ cos(theta_c2_right); sin(theta_c2_right) ]';

xy_inc1_right0 = zeros(n_inc1*ntheta,2);
xy_inc2_right0 = zeros(n_inc2*ntheta,2);
for j = 1:n_inc1
    xy_inc1_right0((j-1)*ntheta+1:j*ntheta,:) = [ xy_c1_right(j,1) + R1 * cos(theta);  xy_c1_right(j,2) + R1 * sin(theta) ]';
    xy_inc2_right0((j-1)*ntheta+1:j*ntheta,:) = [ xy_c2_right(j,1) + R2 * cos(theta);  xy_c2_right(j,2) + R2 * sin(theta) ]';
end

xy_inc_right0 = [ xy_inc1_right0; xy_inc2_right0 ];
len = size(xy_inc_right0,1);

% all the inclusions in the left super cells
xy_inc_right = zeros(n_sup*len,2);
for n = 1:n_sup
    xy_inc_right((n-1)*len+1:n*len,:) = [(n-1)*lat_c, 0] + xy_inc_right0;
end
