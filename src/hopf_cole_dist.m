function phi = hopf_cole_dist(F, V, dirichlet_E, neumann_E, time_step)
%average_edge = compute_average_edge_length(F, V);
%dt = average_edge * average_edge;
%dt = dt * time_step;
dt = time_step;

% solve screened poisson problem
u = solve_screened_poisson(F, V, dirichlet_E, neumann_E, dt);
v = 1-u;

% d(x) = lim t->0 -sqrt(t) ln(u(x))
phi = -sqrt(dt) .* real(log(v));
end


function c = compute_threshold(E, phi)
% compute c = mean(phi(pi)) where pi are points on the input surface
vi = E(:);
vi = unique(vi(:));
c = mean(phi(vi(:)));
end


function avg_edge = compute_average_edge_length(F, V)
element_edge_lengths = compute_edge_length(F, V);
avg_edge = mean(element_edge_lengths);
end


function edge_lengths = compute_edge_length(F, V)
edges = zeros(size(F,1), 3*2);
% all edges of one F
edges(:, 1) = F(:, 1);
edges(:, 2) = F(:, 2);
edges(:, 3) = F(:, 2);
edges(:, 4) = F(:, 3);
edges(:, 5) = F(:, 3);
edges(:, 6) = F(:, 1);

edges_length3 = zeros(size(F,1), 3);
edges_length3(:, 1) = compute_length(V(edges(:,1),:), V(edges(:,2),:));
edges_length3(:, 2) = compute_length(V(edges(:,3),:), V(edges(:,4),:));
edges_length3(:, 3) = compute_length(V(edges(:,5),:), V(edges(:,6),:));

edge_lengths = edges_length3(:);
end


function l = compute_length(p1, p2)
sqr_len = (p2(:,1) - p1(:,1)).^2 + (p2(:,2) - p1(:,2)).^2;
l = sqrt(sqr_len);
end
