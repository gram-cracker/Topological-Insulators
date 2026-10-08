
function [p,t,t_num] = mesh_generation_C4_2PC_PEC(xy_cell,xy_inc_left,xy_inc_right,n_inc,n_sup)

% generate the mesh

% Input
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


% -------------------------------------------------------------------------
% Define vertices, edges, and subregions
% -------------------------------------------------------------------------


%%%%%%%%%
% Nodes %
%%%%%%%%%

nodes = [ xy_inc_left; xy_inc_right; xy_cell ];
N_cell = size(xy_cell,1);

N_inc_left = size(xy_inc_left,1);
N_inc0_left = N_inc_left/(n_inc*n_sup);

N_inc_right = size(xy_inc_right,1);
N_inc0_right = N_inc_right/(n_inc*n_sup);

%%%%%%%%%%
% Edges  %
%%%%%%%%%%

% edges for the left inclusions
edges_inc_left = zeros(N_inc_left,2);

for n = 1:n_sup
    idx = (n-1)*N_inc0_left*n_inc;
    for j = 1:n_inc
        edges_inc_left(idx+(j-1)*N_inc0_left+1:idx+j*N_inc0_left,1) = idx + (j-1)*N_inc0_left + (1:N_inc0_left);
        edges_inc_left(idx+(j-1)*N_inc0_left+1:idx+j*N_inc0_left,2) = idx + (j-1)*N_inc0_left + [2:N_inc0_left, 1];
    end
end

% edges for the right inclusions
edges_inc_right = zeros(N_inc_right,2);

for n = 1:n_sup
    idx = (n-1)*N_inc0_right*n_inc;
    for j = 1:n_inc
        edges_inc_right(idx+(j-1)*N_inc0_right+1:idx+j*N_inc0_right,1) = N_inc_left + idx + (j-1)*N_inc0_right + (1:N_inc0_right);
        edges_inc_right(idx+(j-1)*N_inc0_right+1:idx+j*N_inc0_right,2) = N_inc_left + idx + (j-1)*N_inc0_right + [2:N_inc0_right, 1];
    end
end

edges_cell = zeros(N_cell,2);
edges_cell(1:N_cell,1) = N_inc_left + N_inc_right + (1:N_cell);
edges_cell(1:N_cell,2) = N_inc_left + N_inc_right + [2:N_cell,1];


% assemble all the edges
edges = [ edges_inc_left; edges_inc_right; edges_cell ];

%%%%%%%%%%%%%%
% subregions %
%%%%%%%%%%%%%%

N_region = 1;
subregions = cell(1,N_region);
subregions(1) = {1:N_inc_left+N_inc_right+N_cell};
       

% -------------------------------------------------------------------------
% Specify mesh size and generate Mesh
% -------------------------------------------------------------------------
hdata.hmax = 0.05;         % maximum meshsize parameter
options.output = true;  % suppress displaying the mesh or not

% generate the mesh
[p,t,t_num] = meshfaces(nodes, edges, subregions, hdata, options);
