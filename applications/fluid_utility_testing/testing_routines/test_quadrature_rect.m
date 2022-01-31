
% Quadrature using quadratic tetrahedra on square plate

close all;
clear all;
clc;

addpath('../../../fluid_function_utils');

% Define integration point locations in area coordinates
int_loc1 = [2/3 1/6 1/6];
int_loc2 = [1/6 2/3 1/6];
int_loc3 = [1/6 1/6 2/3];

side = 0.2;
area = 0.2*0.2;
%fun = @(x,y) ((x./side).^5).*((y./side).^3);
%fun = @(x,y) ((x./side).^4).*exp(y./side);
fun = @(x,y) 1.*(x.^0).*(y.^0);

elem_h_list = [];
integrals = [];
for l = 1:1
  % Import structural data
  if (l==1)
    elem_node = importdata('../fem_data/new_circ_data/rect_plate/plate_elems.txt'); elem_node = elem_node.data;
    node_elem = dlmread('../fem_data/new_circ_data/rect_plate/plate_node_elem.txt');
    nodes = dlmread('../fem_data/new_circ_data/rect_plate/plate_nodes.txt');
    structural_node_locs = nodes(:,2:4);
    surface_node_data = dlmread('../fem_data/new_circ_data/rect_plate/plate_face_nodes.txt');
    struct_faces = dlmread('../fem_data/new_circ_data/rect_plate/surf_faces.txt');
  elseif (l==2)
    elem_node = importdata('../fem_data/new_circ_data/rect_plate_mid/rect_mid_elements.txt'); elem_node = elem_node.data;
    node_elem = dlmread('../fem_data/new_circ_data/rect_plate_mid/rect_mid_node_elem.txt');
    nodes = dlmread('../fem_data/new_circ_data/rect_plate_mid/rect_mid_nodes.txt');
    structural_node_locs = nodes(:,2:4);
    surface_node_data = dlmread('../fem_data/new_circ_data/rect_plate_mid/rect_mid_face_nodes.txt');
    struct_faces = dlmread('../fem_data/new_circ_data/rect_plate_mid/rect_mid_faces.txt');
  elseif (l==3)
    elem_node = importdata('../fem_data/new_circ_data/rect_plate_fine/rect_fine_elements.txt'); elem_node = elem_node.data;
    node_elem = dlmread('../fem_data/new_circ_data/rect_plate_fine/rect_fine_node_elem.txt');
    nodes = dlmread('../fem_data/new_circ_data/rect_plate_fine/rect_fine_nodes.txt');
    structural_node_locs = nodes(:,2:4);
    surface_node_data = dlmread('../fem_data/new_circ_data/rect_plate_fine/rect_fine_face_nodes.txt');
    struct_faces = dlmread('../fem_data/new_circ_data/rect_plate_fine/rect_fine_faces.txt');
  end

  elem_a = area / max(size(struct_faces));
  elem_h_list = [elem_h_list, sqrt(elem_a*4/sqrt(3))];

  integral_val = 0;

  counter = zeros(1,4);

  N = max(size(struct_faces));
  for i = 1:1

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

    disp(A1)
    local_integration = quadrature(integrand1*A1, integrand2*A2, integrand3*A3);
    integral_val = integral_val + local_integration;

  end

  integrals = [integrals, integral_val];

end

integral_valid = integral2(fun, 0, 0.2, 0, 0.2);

diffs = abs(integrals - integral_valid);
%c = 10^(log10(diffs(end)) - 3*log10(elem_h_list(end)));
c = diffs(end) / (elem_h_list(end)^3);
yc = c.*elem_h_list.^3;
figure;
loglog(elem_h_list, diffs, '-ob'); hold on;
loglog(elem_h_list, yc, '--r'); grid on;
legend('Numerical', '$\mathcal{O}(h^3)$', 'interpreter', 'latex', 'fontsize', 18);
xlabel('$h$', 'interpreter', 'latex', 'fontsize', 18);
ylabel('$|| I_n - I_a ||$', 'interpreter', 'latex', 'fontsize', 18);
title('Numerical Integration Error: $(x/L)^4 \exp (y/L)$', 'interpreter', 'latex', 'fontsize', 18);


function [int] = quadrature(val1, val2, val3)

  w1 = 1/6;
  w2 = 1/6;
  w3 = 1/6;

  int = w1*val1 + w2*val2 + w3*val3;

end

