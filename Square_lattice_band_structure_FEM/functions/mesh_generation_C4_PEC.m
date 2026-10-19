
function [p,t,t_num,idx] = mesh_generation_C4_PEC(xy_cell,xy_inc)

% generate the mesh

% Input
%
%   xy_cell: grid points for the periodic cell
%
%   xy_inc:grid points for the inclusion
%
%   xy_box: the whole rectangular computational domain
%
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

nodes = [ xy_inc; xy_cell];
N_cell = size(xy_cell,1);
N_inc = size(xy_inc,1);

%%%%%%%%%%
% Edges  %
%%%%%%%%%%

% edges for the inclusion
edges_inc = zeros(N_inc,2);
edges_inc(1:N_inc,1) = 1:N_inc;
edges_inc(1:N_inc,2) = [2:N_inc, 1];

% edges for the cell
edges_cell = zeros(N_cell,2);
edges_cell(1:N_cell,1) = N_inc+(1:N_cell);
edges_cell(1:N_cell,2) = N_inc+[2:N_cell,1];


% assemble all the edges
edges = [ edges_inc; edges_cell ];

%%%%%%%%%%%%%%
% subregions %
%%%%%%%%%%%%%%

N_region = 1;
subregions = cell(1,N_region);

subregions(1) = {1:N_inc+N_cell};
       

% -------------------------------------------------------------------------
% Specify mesh size and generate Mesh
% -------------------------------------------------------------------------
hdata.hmax = 0.02;         % maximum meshsize parameter
options.output = true;  % suppress displaying the mesh or not

% generate the mesh
[p,t,t_num] = meshfaces(nodes, edges, subregions, hdata, options);
