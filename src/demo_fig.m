clear;
clc;
close all;


in_basename = '../data/2d/disk';
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

figure;
title('Hopf-Cole-Varadhan');
trisurf(F, V(:,1), V(:,2), zeros(size(V,1),1), 'CData', phi(:), 'EdgeColor', 'none');
xc = 0.5*(min(V(:,1))+max(V(:,1)));
yc = 0.5*(min(V(:,2))+max(V(:,2)));
campos([xc,yc,1]);
camtarget([xc,yc,0]);
axis equal;
colormap('jet');
shading interp;

figure; 
title('Hopf-Cole-Varadhan');
trisurf(F, V(:,1), V(:,2), phi(:));
colormap('jet');
axis equal;
shading interp;
