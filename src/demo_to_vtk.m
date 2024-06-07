clear;
clc;
close all;


in_basename = '../data/2d/disk';
out_basename = 'disk';
time_step = 0.1;


% Read input data 
ele_name = strcat(in_basename, '.1.ele');
F = read_elements(ele_name);

node_name = strcat(in_basename, '.1.node');
V = read_nodes(node_name);

poly_name = strcat(in_basename, '.1.poly');
E = read_poly(poly_name);

% Boundary conditions:
dirichlet_E = E;
neumann_E = [];

% Compute the dist by solving a screened Poisson problem and 
% using the Hopf-Cole transform / Varahdan's formula
t = F;
p = V;
% Use the Varadhan approach
phi = hopf_cole_dist(t, p(:,1:2), dirichlet_E, neumann_E, time_step);

% Save result as a 2d field and a height field using VTK file format
field_name = strcat(out_basename, '.vtk');
save_2d_field_as_vtk(V, F, phi, 'node', field_name, 'Hopf-Cole-Varadhan');

height_name = strcat(out_basename, '_heightfield.vtk');
save_heightfield_as_vtk(V, F, phi, height_name, 'Hopf-Cole-Varadhan');
