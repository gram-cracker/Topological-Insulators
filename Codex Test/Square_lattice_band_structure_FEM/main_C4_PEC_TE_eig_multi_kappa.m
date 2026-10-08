
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                     %
%   eigensolver for square lattice photoic crystal with C4 symmetry   %
%   Dirichelt boundary condition is imposed on the inclusion boundary %
%   multiple Bloch vector                                             %
%                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

addpath('./mesh2d');
addpath('./finite_elements');
addpath('./functions');
warning('off','all')


% -------------------------------------------------------------------------
% Geometry description:
% -------------------------------------------------------------------------

% periodic cell

lat_c = 1;     % square lattice constant
Lx = lat_c;     Ly = lat_c;
hx = 0.02;   hy = 0.02;      % (smaller than the mesh size hmax)
Nx = Lx/hx; Ny = Ly/hy;
xy_cell = [ zeros(1,Ny), hx:hx:Lx,    Lx*ones(1,Ny), Lx-hx:-hx:0; ...
            Ly-hy:-hy:0, zeros(1,Nx), hy:hy:Ly,      Ly*ones(1,Nx), ]';

xy_box = [ 0, Lx, 0, Ly ];
xy_center = [ Lx/2; Ly/2 ];


% geometry of the inclusions
geo_PEC_inclusions

% -------------------------------------------------------------------------
% generate the mesh
% -------------------------------------------------------------------------
[p,t,t_num] = mesh_generation_C4_PEC(xy_cell,xy_inc);

% -------------------------------------------------------------------------
% build FEM matrics
% -------------------------------------------------------------------------
[M,Mx,My,Me,idx] = FEM_matrices_C4_PEC_TE(p,t,xy_box,xy_inc);


% -------------------------------------------------------------------------
% 1. Solve the eigenvalues for the Bloch wave vectors along the M->X->T->M
% -------------------------------------------------------------------------
% Bloch wave vectors
n_kappa = 100;
%kappa_vec1 = 1/lat_c * [ linspace(0,pi,n_kappa)',  zeros(n_kappa,1) ];         % T -> X
%kappa_vec2 = 1/lat_c * [ pi * ones(n_kappa,1),   linspace(0,pi,n_kappa)' ];    % X -> M
%kappa_vec3 = 1/lat_c * [ linspace(pi,0, n_kappa)', linspace(pi,0, n_kappa)' ]; % M -> T

%kappa_vec = [ kappa_vec1; kappa_vec2; kappa_vec3 ];
% kappa_vec = kappa_vec3;


% kappa_vec1 = 1/lat_c * [ pi * ones(n_kappa,1),   linspace(pi,0,n_kappa)' ];  % M -> X
% kappa_vec2 = 1/lat_c * [ linspace(pi,0,n_kappa)',  zeros(n_kappa,1) ];       % X -> T
% kappa_vec3 = 1/lat_c * [ linspace(0,pi,n_kappa)', linspace(0,pi,n_kappa)' ]; % T -> M
%
 %kappa_vec = [ kappa_vec1; kappa_vec2; kappa_vec3 ];

 kappa_vec1 = [ pi * ones(n_kappa,1),   linspace(pi,0,n_kappa)' ];  % X -> N
 kappa_vec2 = [ linspace(pi,0,n_kappa)',  zeros(n_kappa,1) ];       % N -> T
 kappa_vec3 = [ zeros(n_kappa,1),   linspace(0,pi,n_kappa)' ];      % T -> M
 kappa_vec4 = [ linspace(0,pi,n_kappa)', pi * ones(n_kappa,1) ];    % M -> X
%
 kappa_vec = 1/lat_c * [ kappa_vec1; kappa_vec2; kappa_vec3; kappa_vec4 ];


% -------------------------------------------------------------------------
% solve the first n_omega eigenvalues over each Bloch wave vector
% -------------------------------------------------------------------------
n_omega = 7;
omega = zeros(n_omega,size(kappa_vec,1));
for n = 1:size(kappa_vec,1)
    kappa = kappa_vec(n,:);
    A = M - 2i * ( kappa(1)*Mx + kappa(2)*My ) + (kappa*kappa')*Me;
    B = Me;
    lambda = eigs(A,B,n_omega,'smallestabs');
    omega(:,n) = sqrt(lambda)*lat_c/(2*pi);
end

figure;
for n = 1:n_omega
    hold on; plot(real(omega(n,:)));
end



% -------------------------------------------------------------------------
% 2. Solve the eigenvalues for Bloch wave vectors over the rectangular grids
% -------------------------------------------------------------------------


n_kappa = 21;
kappa = [0,0];
kappa_vec1 = kappa(1) + 0.1*pi/lat_c * linspace(-1,1,n_kappa);
kappa_vec2 = kappa(2) + 0.1*pi/lat_c * linspace(-1,1,n_kappa);
[ Kappa1, Kappa2 ] = meshgrid(kappa_vec1,kappa_vec2);

n_omega = 5;
omega = zeros(n_kappa,n_kappa,n_omega);

for m = 1:n_kappa
    m
    for n = 1:n_kappa
        kappa = [ Kappa1(m,n), Kappa2(m,n) ] ;
        A = M - 2i * ( kappa(1)*Mx + kappa(2)*My ) + (kappa*kappa')*Me;
        B = Me;
        lambda = eigs(A,B,n_omega,'smallestabs');
        omega(m,n,:) = sqrt(lambda); % lat_c*sqrt(lambda)/(2*pi);
    end
end

figure;
for n = 1:n_omega
    hold on; surf(Kappa1/pi,Kappa2/pi,real(omega(:,:,n)));
end
