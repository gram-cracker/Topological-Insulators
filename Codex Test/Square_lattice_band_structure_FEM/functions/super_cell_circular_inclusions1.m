

% ---------------------------------
% circular inclusions in the cells
% ---------------------------------

n_inc = 4;                                      % number of inclusions
theta_c = pi/4 + 2*pi/n_inc * ( 0:1:n_inc-1 ); 

ntheta_left = 20;
dtheta_left = 2*pi/ntheta_left;
theta_left = 0:dtheta_left:2*pi-dtheta_left;

% inclusion in the first left cell
R_left = 0.14 * lat_c;                       % radius of the inclusions
d_left = sqrt(2) * 2.5 * R_left;             % distance between the orgin and the inclusion center
xyc_left = repmat([ lat_c/2, lat_c/2 ],4,1) + d_left * [ cos(theta_c); sin(theta_c) ]';

xy_inc0_left  = zeros(n_inc*ntheta_left,2);
for j = 1:n_inc
    xy_inc0_left((j-1)*ntheta_left+1:j*ntheta_left,:) = [ xyc_left(j,1) + R_left * cos(theta_left);  xyc_left(j,2) + R_left * sin(theta_left) ]';
end

% all the inclusions in the left super cells
xy_inc_left = zeros(4*n_sup*ntheta_left,2);
for n = 1:n_sup
    xy_inc_left((n-1)*4*ntheta_left+1:n*4*ntheta_left,:) = [-n*lat_c, 0] + xy_inc0_left;
    % hold on; plot(xy_inc_left((n-1)*4*ntheta_left+1:(4*n-3)*ntheta_left,1), xy_inc_left((n-1)*4*ntheta_left+1:(4*n-3)*ntheta_left,2) );
    % hold on; plot(xy_inc_left((4*n-3)*ntheta_left+1:(4*n-2)*ntheta_left,1), xy_inc_left((4*n-3)*ntheta_left+1:(4*n-2)*ntheta_left,2) );
    % hold on; plot(xy_inc_left((4*n-2)*ntheta_left+1:(4*n-1)*ntheta_left,1), xy_inc_left((4*n-2)*ntheta_left+1:(4*n-1)*ntheta_left,2) );
    % hold on; plot(xy_inc_left((4*n-1)*ntheta_left+1:(4*n)*ntheta_left,1), xy_inc_left((4*n-1)*ntheta_left+1:(4*n)*ntheta_left,2) );
end


ntheta_right = 20;
dtheta_right = 2*pi/ntheta_right;
theta_right = 0:dtheta_right:2*pi-dtheta_right;

% inclusion in the first right cell
R_right = 0.06 * lat_c;                        % radius of the inclusions
d_right = sqrt(2) * 2.5 * R_right;             % distance between the orgin and the inclusion center
xyc_right = repmat([ lat_c/2, lat_c/2 ],4,1) + d_right * [ cos(theta_c); sin(theta_c) ]';

xy_inc0_right  = zeros(n_inc*ntheta_right,2);
for j = 1:n_inc
    xy_inc0_right((j-1)*ntheta_right+1:j*ntheta_right,:) = [ xyc_right(j,1) + R_right * cos(theta_right);  xyc_right(j,2) + R_right * sin(theta_right) ]';
end

% all the inclusions in the right super cells
xy_inc_right = zeros(4*n_sup*ntheta_right,2);
for n = 1:n_sup
    xy_inc_right((n-1)*4*ntheta_right+1:n*4*ntheta_right,:) = [(n-1)*lat_c, 0] + xy_inc0_right;
    % hold on; plot(xy_inc_right((n-1)*4*ntheta_right+1:(4*n-3)*ntheta_right,1), xy_inc_right((n-1)*4*ntheta_right+1:(4*n-3)*ntheta_right,2) );
    % hold on; plot(xy_inc_right((4*n-3)*ntheta_right+1:(4*n-2)*ntheta_right,1), xy_inc_right((4*n-3)*ntheta_right+1:(4*n-2)*ntheta_right,2) );
    % hold on; plot(xy_inc_right((4*n-2)*ntheta_right+1:(4*n-1)*ntheta_right,1), xy_inc_right((4*n-2)*ntheta_right+1:(4*n-1)*ntheta_right,2) );
    % hold on; plot(xy_inc_right((4*n-1)*ntheta_right+1:(4*n)*ntheta_right,1), xy_inc_right((4*n-1)*ntheta_right+1:(4*n)*ntheta_right,2) );
end
