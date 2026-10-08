
function [p,t,t_num,neu_nodes] = mesh_generation_C4_scat(xy_box,xy_cell,xy_inc_bot,n_inc,n_sup)

% generate the mesh

% Input
%
%   xy_box: double, [x_min, x_max, y_min, y_max] denote the domain
%
%   xy_cell: grid points for the periodic cell
%
%   xy_inc: grid points for the inclusions in the supercell
% 
%   n_inc: number of inclusions in each cell
%
%   n_sup: number of supercells
%
% Output
%
%   p: double, (2,n_nodes) array of mesh nodes
%
%   t: int, (n_elements, n_dof) triangular element matrix
%
%   t_num: array of regions numbers for each triangle in t
%
%   neu_nodes: Neumann nodes 


x_min = xy_box(1); x_max = xy_box(2); 
y_min = xy_box(3); y_max = xy_box(4);


% -------------------------------------------------------------------------
% Define vertices, edges, and subregions
% -------------------------------------------------------------------------


%%%%%%%%%
% Nodes %
%%%%%%%%%

nodes = [ xy_inc_bot; xy_cell ];
N_cell = size(xy_cell,1);

N_inc_bot = size(xy_inc_bot,1);
N_inc0_bot = N_inc_bot/(n_inc*n_sup);

%%%%%%%%%%
% Edges  %
%%%%%%%%%%

% edges for the bottom inclusions
edges_inc_bot = zeros(N_inc_bot,2);

for n = 1:n_sup
    idx = (n-1)*N_inc0_bot*n_inc;
    for j = 1:n_inc
        edges_inc_bot(idx+(j-1)*N_inc0_bot+1:idx+j*N_inc0_bot,1) = idx + (j-1)*N_inc0_bot + (1:N_inc0_bot);
        edges_inc_bot(idx+(j-1)*N_inc0_bot+1:idx+j*N_inc0_bot,2) = idx + (j-1)*N_inc0_bot + [2:N_inc0_bot, 1];
    end
end

edges_cell = zeros(N_cell,2);
edges_cell(1:N_cell,1) = N_inc_bot + (1:N_cell);
edges_cell(1:N_cell,2) = N_inc_bot + [2:N_cell,1];


% assemble all the edges
edges = [ edges_inc_bot; edges_cell ];

%%%%%%%%%%%%%%
% subregions %
%%%%%%%%%%%%%%

N_region = 2;
subregions = cell(1,N_region);

subregions(1) = {1:N_inc_bot};
subregions(2) = {1:N_inc_bot+N_cell};
       

% -------------------------------------------------------------------------
% Specify mesh size and generate Mesh
% -------------------------------------------------------------------------
hdata.hmax = 0.025;      % maximum meshsize parameter
options.output = true;  % suppress displaying the mesh or not

% generate the mesh
[p,t,t_num] = meshfaces(nodes, edges, subregions, hdata, options);


% -------------------------------------------------------------------------
% Locate Neumann nodes (upper boundary)
% -------------------------------------------------------------------------
tol = 1e-8*hdata.hmax; 
n_nodes = size(p,1);
i_all = (1:n_nodes)';
i_neu = i_all(abs(p(:,2)-y_max)<tol);  % nodes on upper bnd
[~,i_srt] = sort(p(i_neu,1));          % sort in x-direction
neu_nodes = i_neu(i_srt);              % neumann nodes, sorted.