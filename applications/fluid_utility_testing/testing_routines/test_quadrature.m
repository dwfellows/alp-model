
close all;
clear all;
clc;

addpath('../../../fluid_function_utils');

% Define integration point locations in area coordinates
int_loc1 = [2/3 1/6 1/6];
int_loc2 = [1/6 2/3 1/6];
int_loc3 = [1/6 1/6 2/3];

R = 0.1;
area = pi*R^2;

fun_num = @(r,t) ((r./R).^3) .* cos(t./(2*pi)).^2;
fun = @(r,t) ( ((r./R).^3) .* (cos(t ./(2*pi)).^2) ) .*r.*(t.^0);

elem_h_list = [];
integrals = [];
for l = 1:5
  % Import structural data
  if (l==1)
    elem_node = dlmread('../fem_data/circular_coarse/coarse.csv');
    node_elem = dlmread('../fem_data/circular_coarse/coarse_node_elem.csv');
    xdir_deform = dlmread('../fem_data/circular_coarse/coarse_xdir_deformation.txt');
    ydir_deform = dlmread('../fem_data/circular_coarse/coarse_ydir_deformation.txt');
    zdir_deform = dlmread('../fem_data/circular_coarse/coarse_zdir_deformation.txt');
    def_shape_function = [xdir_deform(:,5), ydir_deform(:,5), zdir_deform(:,5)];
    def_shape_function_norms = vecnorm(def_shape_function'); max_def = max(def_shape_function_norms);
    def_shape_function = def_shape_function ./ max_def;
    structural_node_locs = xdir_deform(:,2:4);
    surface_node_data = dlmread('../fem_data/circular_coarse/coarse_surface.txt');
    struct_faces = dlmread('../fem_data/circular_coarse/struc_data/coarse_faces.txt');
  elseif (l==3)
    elem_node = dlmread('../fem_data/circular_mid/mid.csv');
    node_elem = dlmread('../fem_data/circular_mid/mid_node_elem.csv');
    xdir_deform = dlmread('../fem_data/circular_mid/mid_xdir_deformation.txt');
    ydir_deform = dlmread('../fem_data/circular_mid/mid_ydir_deformation.txt');
    zdir_deform = dlmread('../fem_data/circular_mid/mid_zdir_deformation.txt');
    def_shape_function = [xdir_deform(:,5), ydir_deform(:,5), zdir_deform(:,5)];
    def_shape_function_norms = vecnorm(def_shape_function'); max_def = max(def_shape_function_norms);
    def_shape_function = def_shape_function ./ max_def;
    structural_node_locs = xdir_deform(:,2:4);
    surface_node_data = dlmread('../fem_data/circular_mid/mid_surface.txt');
    struct_faces = dlmread('../fem_data/circular_mid/struc_data/mid_faces.txt');
  elseif (l==5)
    elem_node = dlmread('../fem_data/circular_fine/fine.csv');
    node_elem = dlmread('../fem_data/circular_fine/fine_node_elem.csv');
    xdir_deform = dlmread('../fem_data/circular_fine/fine_xdir_deformation.txt');
    ydir_deform = dlmread('../fem_data/circular_fine/fine_ydir_deformation.txt');
    zdir_deform = dlmread('../fem_data/circular_fine/fine_zdir_deformation.txt');
    def_shape_function = [xdir_deform(:,5), ydir_deform(:,5), zdir_deform(:,5)];
    def_shape_function_norms = vecnorm(def_shape_function'); max_def = max(def_shape_function_norms);
    def_shape_function = def_shape_function ./ max_def;
    structural_node_locs = xdir_deform(:,2:4);
    surface_node_data = dlmread('../fem_data/circular_fine/fine_surface.txt');
    struct_faces = dlmread('../fem_data/circular_fine/struc_data/fine_faces.txt');
  elseif (l==2)
    %elem_node = dlmread('../fem_data/circular_cm/cm_elements.txt');
    elem_node = importdata('../fem_data/circular_cm/cm_elements.txt');
    elem_node = elem_node.data;
    node_elem = dlmread('../fem_data/circular_cm/cm_node_elem.csv');
    xdir_deform = dlmread('../fem_data/circular_cm/cm_xdir_deformation.txt');
    ydir_deform = dlmread('../fem_data/circular_cm/cm_ydir_deformation.txt');
    zdir_deform = dlmread('../fem_data/circular_cm/cm_zdir_deformation.txt');
    def_shape_function = [xdir_deform(:,5), ydir_deform(:,5), zdir_deform(:,5)];
    def_shape_function_norms = vecnorm(def_shape_function'); max_def = max(def_shape_function_norms);
    def_shape_function = def_shape_function ./ max_def;
    structural_node_locs = xdir_deform(:,2:4);
    surface_node_data = dlmread('../fem_data/circular_cm/cm_surf_nodes.txt');
    struct_faces = dlmread('../fem_data/circular_cm/struc_data/cm_faces.txt');
  elseif (l==4)
    %elem_node = dlmread('../fem_data/circular_fm/fm_elements.txt');
    elem_node = importdata('../fem_data/circular_fm/fm_elements.txt');
    elem_node = elem_node.data;
    node_elem = dlmread('../fem_data/circular_fm/fm_node_elem.csv');
    xdir_deform = dlmread('../fem_data/circular_fm/fm_xdir_deformation.txt');
    ydir_deform = dlmread('../fem_data/circular_fm/fm_ydir_deformation.txt');
    zdir_deform = dlmread('../fem_data/circular_fm/fm_zdir_deformation.txt');
    def_shape_function = [xdir_deform(:,5), ydir_deform(:,5), zdir_deform(:,5)];
    def_shape_function_norms = vecnorm(def_shape_function'); max_def = max(def_shape_function_norms);
    def_shape_function = def_shape_function ./ max_def;
    structural_node_locs = xdir_deform(:,2:4);
    surface_node_data = dlmread('../fem_data/circular_fm/fm_surface.txt');
    struct_faces = dlmread('../fem_data/circular_fm/struc_data/fm_faces.txt');
  end

  elem_a = area / max(size(struct_faces));
  elem_h_list = [elem_h_list, sqrt(elem_a*4/sqrt(3))];

  integral_val = 0;

  N = max(size(struct_faces));
  for i = 1:N

    % Acquire local element, nodes, surface
    loc_elem = struct_faces(i,1);
    loc_elem_surf_nodes = struct_faces(i,2:11);

    loc_elem_nodes = elem_node(loc_elem,:);
    loc_elem_nodelocs = structural_node_locs(loc_elem_nodes,:);
    loc_elem_surf_num = find(loc_elem_surf_nodes(1:4) == 0);
    loc_elem_surf = [1 2 3 4]; loc_elem_surf(loc_elem_surf_num) = [];

    % Define surface quantity on nodes (analytical)
    int_nodes = zeros(10,1);
    for m = 1:10
      x = loc_elem_nodelocs(m,1); y = loc_elem_nodelocs(m,2);
      r = sqrt(x^2 + y^2); t = atan2(y,x); t = wrapTo2Pi(t);
      int_nodes(m) = fun_num(r,t);
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

integral_valid = integral2(fun, 0, R, 0, 2*pi);

diffs = abs(integrals - integral_valid);
%c = 10^(log10(diffs(end)) - 3*log10(elem_h_list(end)));
c = diffs(end) / (elem_h_list(end)^3);
yc = c.*elem_h_list.^3;
figure;
loglog(elem_h_list, diffs, '-ob'); hold on;
loglog(elem_h_list, yc, '--r'); grid on;
legend('Numerical', '$\mathcal{O}(h^3)$', 'interpreter', 'latex', 'location', 'southwest', 'fontsize', 18);
xlabel('$h$', 'interpreter', 'latex', 'fontsize', 18);
ylabel('$| I_n - I_a |$', 'interpreter', 'latex', 'fontsize', 18);
title('Numerical Integration Error: $f(r) = (r/R)^3 \cdot \cos^2 (t / (2\pi))$', 'interpreter', 'latex', 'fontsize', 18);


function [int] = quad(val1, val2, val3)

  w1 = 1/6;
  w2 = 1/6;
  w3 = 1/6;

  int = w1*val1 + w2*val2 + w3*val3;

end


function [returned_nodes] = correct_nodes(given_nodes)

    returned_nodes(1:4,:) = given_nodes(1:4,:);

    returned_nodes(5,:) = (given_nodes(1,:) + given_nodes(2,:)) ./ 2;
    returned_nodes(6,:) = (given_nodes(2,:) + given_nodes(3,:)) ./ 2;
    returned_nodes(7,:) = (given_nodes(3,:) + given_nodes(1,:)) ./ 2;
    returned_nodes(8,:) = (given_nodes(1,:) + given_nodes(4,:)) ./ 2;
    returned_nodes(9,:) = (given_nodes(2,:) + given_nodes(4,:)) ./ 2;
    returned_nodes(10,:) = (given_nodes(3,:) + given_nodes(4,:)) ./ 2;


end
