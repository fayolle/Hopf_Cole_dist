function u = solve_screened_poisson(F, V, dirichlet, neumann, time_step)
% preparation of data
free_nodes = setdiff(1:size(V,1),unique(dirichlet));
L = sparse(size(V,1),size(V,1));
M = sparse(size(V,1),size(V,1));

dt = time_step;

% cotan weight matrix
for j = 1:size(F,1)
    L(F(j,:),F(j,:)) = L(F(j,:), ...
        F(j,:)) + compute_local_assembly_matrix(V(F(j,:),:));
end

% mass matrix
for j = 1:size(F,1)
    M(F(j,:),F(j,:)) = M(F(j,:), ...
        F(j,:)) + det([1,1,1;V(F(j,:),:)'])...
        *[2,1,1;1,2,1;1,1,2]/24;
end

b = sparse(size(V,1),1);

% Volume Forces
for j = 1:size(F,1)
    b(F(j,:)) = b(F(j,:)) + ...
        det([1,1,1; V(F(j,:),:)']) * ...
        compute_volume_force(sum(V(F(j,:),:))/3)/6;
end

% Neumann conditions
for j = 1 : size(neumann,1)
    b(neumann(j,:)) = b(neumann(j,:)) + ...
        norm(V(neumann(j,1),:)-V(neumann(j,2),:))*...
        compute_neumann(sum(V(neumann(j,:),:))/2)/2;
end

% Dirichlet conditions
u = sparse(size(V,1),1);
u(unique(dirichlet)) = compute_dirichlet(V(unique(dirichlet),:));
b = b - (dt * L + M) * u;

% Computation of the solution
u(free_nodes) = (dt*L(free_nodes, free_nodes)+ ...
    M(free_nodes, free_nodes))\b(free_nodes);

end


function L = compute_local_assembly_matrix(vertices)
d = size(vertices,2);
G = [ones(1,d+1);vertices'] \ [zeros(1,d);eye(d)];
L = det([ones(1,d+1);vertices']) * (G*G') / prod(1:d);
end


function volume_force = compute_volume_force(x)
volume_force = ones(size(x,1),1);
end


function neumann_bc = compute_neumann(x)
neumann_bc = zeros(size(x,1),1);
end


function dirichlet_bc = compute_dirichlet(x)
dirichlet_bc = zeros(size(x,1),1);
end
