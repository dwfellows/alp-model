
% Quadrature using quadratic tetrahedra on sphere surface

close all;
clear all;
clc;

addpath('../../../fluid_function_utils');

fs = 24;

% Define integration point locations in area coordinates
int_loc1 = [2/3 1/6 1/6];
int_loc2 = [1/6 2/3 1/6];
int_loc3 = [1/6 1/6 2/3];

R = 0.2;
area = 4*pi*(R^2);
fun = @(x,y) 1;

elem_h_list = [];
integrals = [];
for l = 1:4
  % Import structural data
  if (l==1)
    elem_node = importdata('../fem_data/sphere/sphere_c1/sphere_c1_elems.txt'); elem_node = elem_node.data;
    node_elem = dlmread('../fem_data/sphere/sphere_c1/sphere_c1_node_elem.txt');
    nodes = dlmread('../fem_data/sphere/sphere_c1/sphere_c1_nodes.txt');
    structural_node_locs = nodes(:,2:4);
    surface_node_data = dlmread('../fem_data/sphere/sphere_c1/sphere_c1_face_nodes.txt');
    struct_faces = dlmread('../fem_data/sphere/sphere_c1/sphere_c1_faces.txt');
  elseif (l==2)
    elem_node = importdata('../fem_data/sphere/sphere_c2/sphere_c2_elems.txt'); elem_node = elem_node.data;
    node_elem = dlmread('../fem_data/sphere/sphere_c2/sphere_c2_node_elem.txt');
    nodes = dlmread('../fem_data/sphere/sphere_c2/sphere_c2_nodes.txt');
    structural_node_locs = nodes(:,2:4);
    surface_node_data = dlmread('../fem_data/sphere/sphere_c2/sphere_c2_face_nodes.txt');
    struct_faces = dlmread('../fem_data/sphere/sphere_c2/sphere_c2_faces.txt');
  elseif (l==3)
    elem_node = importdata('../fem_data/sphere/sphere_c3/sphere_c3_elems.txt'); elem_node = elem_node.data;
    node_elem = dlmread('../fem_data/sphere/sphere_c3/sphere_c3_node_elem.txt');
    nodes = dlmread('../fem_data/sphere/sphere_c3/sphere_c3_nodes.txt');
    structural_node_locs = nodes(:,2:4);
    surface_node_data = dlmread('../fem_data/sphere/sphere_c3/sphere_c3_face_nodes.txt');
    struct_faces = dlmread('../fem_data/sphere/sphere_c3/sphere_c3_faces.txt');
  elseif (l==4)
    elem_node = importdata('../fem_data/sphere/sphere_c4/sphere_c4_elems.txt'); elem_node = elem_node.data;
    node_elem = dlmread('../fem_data/sphere/sphere_c4/sphere_c4_node_elem.txt');
    nodes = dlmread('../fem_data/sphere/sphere_c4/sphere_c4_nodes.txt');
    structural_node_locs = nodes(:,2:4);
    surface_node_data = dlmread('../fem_data/sphere/sphere_c4/sphere_c4_face_nodes.txt');
    struct_faces = dlmread('../fem_data/sphere/sphere_c4/sphere_c4_faces.txt');
  end

  elem_a = area / max(size(struct_faces));
  elem_h_list = [elem_h_list, sqrt(elem_a*4/sqrt(3))];

  integral_val = 0;

  counter = zeros(1,4);

  N = max(size(struct_faces));
  for i = 1:N

    % Acquire local element, nodes, surface
    loc_elem = struct_faces(i,1);
    loc_elem_surf_nodes = struct_faces(i,2:11);

    loc_elem_nodes = elem_node(loc_elem,:);
    loc_elem_nodelocs = structural_node_locs(loc_elem_nodes,:);
    loc_elem_surf_num = find(loc_elem_surf_nodes(1:4) == 0);
    loc_elem_surf = [1 2 3 4]; loc_elem_surf(loc_elem_surf_num) = [];

    counter(loc_elem_surf_num) = counter(loc_elem_surf_num) + 1;

    % Define surface quantity on nodes (analytical)
    int_nodes = zeros(10,1);
    for m = 1:10
      x = loc_elem_nodelocs(m,1); y = loc_elem_nodelocs(m,2);
      int_nodes(m) = fun(x,y);
    end
  
    % Interpolate integrand at integration points
    loc1_bary = zeros(1,4); loc1_bary(loc_elem_surf) = int_loc1;
    loc1 = tet_interp(loc_elem_nodelocs, loc1_bary);
    integrand1 = tet_interp(int_nodes, loc1_bary);
    [~,j1,~,~] = jacobian(loc1_bary(2), loc1_bary(3), loc1_bary(4), loc_elem_nodelocs);
    [dxdLi1, dxdLj1, order1] = def_derivs(loc_elem_nodelocs, int_nodes, loc1, loc1_bary, loc_elem_surf_num);
    A1 = cross(dxdLi1, dxdLj1); A1 = norm(A1);

    tri_j = zeros(2,2);
    tri_j(1,1) = (4*loc1_bary(1)-1)*loc_elem_nodelocs(1,1) + (-3+4*loc1_bary(1)+4*loc1_bary(2))*loc_elem_nodelocs(3,1) ...
               + 4*loc1_bary(2)*loc_elem_nodelocs(5,1) -4*loc1_bary(2)*loc_elem_nodelocs(6,1) + 4*(1-2*loc1_bary(1)-loc1_bary(2))*loc_elem_nodelocs(7,1);


    tri_j(1,2) = (4*loc1_bary(1)-1)*loc_elem_nodelocs(1,2) + (-3+4*loc1_bary(1)+4*loc1_bary(2))*loc_elem_nodelocs(3,2) ...
               + 4*loc1_bary(2)*loc_elem_nodelocs(5,2) -4*loc1_bary(2)*loc_elem_nodelocs(6,2) + 4*(1-2*loc1_bary(1)-loc1_bary(2))*loc_elem_nodelocs(7,2);


    tri_j(2,1) = (4*loc1_bary(2)-1)*loc_elem_nodelocs(2,1) + (-3+4*loc1_bary(1)+4*loc1_bary(2))*loc_elem_nodelocs(3,1) ...
               + 4*loc1_bary(1)*loc_elem_nodelocs(5,1) + 4*(1-loc1_bary(1)-2*loc1_bary(2))*loc_elem_nodelocs(6,1) - 4*loc1_bary(1)*loc_elem_nodelocs(7,1);


    tri_j(2,2) = (4*loc1_bary(2)-1)*loc_elem_nodelocs(2,2) + (-3+4*loc1_bary(1)+4*loc1_bary(2))*loc_elem_nodelocs(3,2) ...
               + 4*loc1_bary(1)*loc_elem_nodelocs(5,2) + 4*(1-loc1_bary(1)-2*loc1_bary(2))*loc_elem_nodelocs(6,2) - 4*loc1_bary(1)*loc_elem_nodelocs(7,2); 


    loc2_bary = zeros(1,4); loc2_bary(loc_elem_surf) = int_loc2;
    loc2 = tet_interp(loc_elem_nodelocs, loc2_bary);
    integrand2 = tet_interp(int_nodes, loc2_bary);
    [~,j2,~,~] = jacobian(loc2_bary(2), loc2_bary(3), loc2_bary(4), loc_elem_nodelocs);
    [dxdLi2, dxdLj2, order2] = def_derivs(loc_elem_nodelocs, int_nodes, loc2, loc2_bary, loc_elem_surf_num);
    A2 = cross(dxdLi2, dxdLj2); A2 = norm(A2);

  
    loc3_bary = zeros(1,4); loc3_bary(loc_elem_surf) = int_loc3;
    loc3 = tet_interp(loc_elem_nodelocs, loc3_bary);
    integrand3 = tet_interp(int_nodes, loc3_bary);
    [~,j3,~,~] = jacobian(loc3_bary(2), loc3_bary(3), loc3_bary(4), loc_elem_nodelocs);
    [dxdLi3, dxdLj3, order3] = def_derivs(loc_elem_nodelocs, int_nodes, loc3, loc3_bary, loc_elem_surf_num);
    A3 = cross(dxdLi3, dxdLj3); A3 = norm(A3);


    local_integration = quad(integrand1*A1, integrand2*A2, integrand3*A3);
    integral_val = integral_val + local_integration;

  end

  integrals = [integrals, integral_val];

end

integral_valid = area;

diffs = abs(integrals - integral_valid) / abs(integral_valid);
%c = 10^(log10(diffs(end)) - 3*log10(elem_h_list(end)));
c = diffs(end) / (elem_h_list(end)^3);
yc = c.*elem_h_list.^3;
figure;
lg1 = loglog(elem_h_list, diffs, '-ob'); hold on;
lg1.LineWidth = 2;
lg2 = loglog(elem_h_list, yc, '--r'); grid on;
lg2.LineWidth = 2;
a = get(gca,'YTickLabel')
set(gca,'YTickLabel',a,'fontsize',18);
legend('Numerical', '$\mathcal{O}(h^3)$', 'interpreter', 'latex', 'location', 'northwest', 'fontsize', 24);
xlabel('Avg. Triangle Side Length [m]', 'interpreter', 'latex', 'fontsize', fs);
ylabel('$|| I_n - I_a || / ||I_a||$', 'interpreter', 'latex', 'fontsize', fs);
title('Numerical Integration Error: Sphere Surface Area', 'interpreter', 'latex', 'fontsize', fs);


function [int] = quad(val1, val2, val3)

  w1 = 1/6;
  w2 = 1/6;
  w3 = 1/6;

  int = w1*val1 + w2*val2 + w3*val3;

end

