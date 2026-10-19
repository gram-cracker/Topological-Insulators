
function [p,t,t_num] = mesh_generation_diamond_PEC(xy_cell,xy_inc1,xy_inc2)

% generate the mesh

% Input
%
%   xy_cell: grid points for the periodic cell
%
%   xy_inc1(2):grid points for the inclusions
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

nodes = [ xy_inc1; xy_inc2; xy_cell];
N_cell = size(xy_cell,1);
N_inc1 = size(xy_inc1,1);
N_inc2 = size(xy_inc2,1);
N_inc = N_inc1 + N_inc2;

%%%%%%%%%%
% Edges  %
%%%%%%%%%%

% edges for the inclusion
edges_inc = zeros(N_inc1+N_inc2,2);
edges_inc(1:N_inc,1) = 1:N_inc;
edges_inc(1:N_inc1,2) = [2:N_inc1, 1];
edges_inc(N_inc1+1:N_inc,2) = [N_inc1+2:N_inc, N_inc1+1];

edges_cell = zeros(N_cell,2);
edges_cell(1:N_cell,1) = N_inc+(1:N_cell);
edges_cell(1:N_cell,2) = N_inc+[2:N_cell,1];

% edges outside the inclusion
% edges_ext = zeros(N_inc+N_cell,2);
% edges_ext(N_inc+(1:N_cell),1) = N_inc+[1:N_cell];
% edges_ext(N_inc+(1:N_cell),2) = N_inc+[2:N_cell, 1];

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
