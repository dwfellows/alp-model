
close all;
clear all;
clc;

fs = 18;

% Add mesh interfacing tools path
addpath('../../../fluid_function_utils');

% Define circular plate values
R = 0.1;
besselfuns = @(l) besselj(0,l)*besseli(1,l) + besseli(0,l)*besselj(1,l);
lambda = fsolve(besselfuns, 3.1962);
k = lambda / R;
h = 0.0001;
growth = 10^0;
area = pi*(R^2);

% Define integration point locations
int_loc1 = [2/3, 1/6, 1/6];
int_loc2 = [1/6, 2/3, 1/6];
int_loc3 = [1/6, 1/6, 2/3];

% Work through random points on all meshes to evaluate convergence
elem_h_list = [];
diff_list = []; loc_diff_list = [];
for l = 1:3

  if (l==1) % rect data
    elem_node = importdata('../fem_data/rect_plate/plate_elems.txt');
    elem_node = elem_node.data;
    node_elem = importdata('../fem_data/rect_plate/plate_node_elem.txt');
    xdir_deform = zeros(max(size(node_elem)), 5);
    ydir_deform = zeros(max(size(node_elem)), 5);
    zdir_deform = zeros(max(size(node_elem)), 5);
    node_locs = dlmread('../fem_data/rect_plate/plate_nodes.txt'); node_locs = node_locs(:,2:4);
    surface_node_data = dlmread('../fem_data/rect_plate/plate_face_nodes.txt');
    struct_faces = dlmread('../fem_data/rect_plate/surf_faces.txt');
    elem_a = area / max(size(struct_faces));
    elem_h_list = [elem_h_list, (elem_a*4/sqrt(3))^(1/2)];
  elseif (l==2) % rect-mid data
    elem_node = importdata('../fem_data/new_circ_data/rect_plate_mid/rect_mid_elements.txt');
    elem_node = elem_node.data;
    node_elem = importdata('../fem_data/new_circ_data/rect_plate_mid/rect_mid_node_elem.txt');
    xdir_deform = zeros(max(size(node_elem)), 5);
    ydir_deform = zeros(max(size(node_elem)), 5);
    zdir_deform = zeros(max(size(node_elem)), 5);
    node_locs = dlmread('../fem_data/new_circ_data/rect_plate_mid/rect_mid_nodes.txt'); node_locs = node_locs(:,2:4);
    surface_node_data = dlmread('../fem_data/new_circ_data/rect_plate_mid/rect_mid_face_nodes.txt');
    struct_faces = dlmread('../fem_data/new_circ_data/rect_plate_mid/rect_mid_faces.txt');
    elem_a = area / max(size(struct_faces));
    elem_h_list = [elem_h_list, (elem_a*4/sqrt(3))^(1/2)];
  elseif (l==3)  % Fine data
    elem_node = importdata('../fem_data/new_circ_data/rect_plate_fine/rect_fine_elements.txt');
    elem_node = elem_node.data;
    node_elem = importdata('../fem_data/new_circ_data/rect_plate_fine/rect_fine_node_elem.txt');
    xdir_deform = zeros(max(size(node_elem)), 5);
    ydir_deform = zeros(max(size(node_elem)), 5);
    zdir_deform = zeros(max(size(node_elem)), 5);
    node_locs = dlmread('../fem_data/new_circ_data/rect_plate_fine/rect_fine_nodes.txt'); node_locs = node_locs(:,2:4);
    surface_node_data = dlmread('../fem_data/new_circ_data/rect_plate_fine/rect_fine_face_nodes.txt');
    struct_faces = dlmread('../fem_data/new_circ_data/rect_plate_fine/rect_fine_faces.txt');
    elem_a = area / max(size(struct_faces));
    elem_h_list = [elem_h_list, (elem_a*4/sqrt(3))^(1/2)];
  end

  max_diff = 0; max_loc_diff = 0;
  counter = [0 0 0 0];
  diffs = [];
  location_diffs = [];


  N = max(size(struct_faces));  
  %% Create following single mode approximation
  % Go over each face, compute numerical normals at each point, compare to analytically known normals
  for i = 1:N

    % Construct integration points
    loc_elem = struct_faces(i,1);
    loc_elem_surf_nodes = struct_faces(i, 2:11);

    loc_elem_nodes = elem_node(loc_elem,:);
    loc_elem_nodelocs = node_locs(loc_elem_nodes,:);
    loc_elem_surf_num = find(loc_elem_surf_nodes(1:4) == 0);

    loc_elem_surf = [1 2 3 4]; loc_elem_surf(loc_elem_surf_num) = [];

    counter(loc_elem_surf_num) = counter(loc_elem_surf_num) + 1;
    
    % Modify local nodes to put mid-nodes on actual mid-line
%    loc_elem_nodelocs(5,:) = (loc_elem_nodelocs(1,:) + loc_elem_nodelocs(2,:)) ./ 2;
%    loc_elem_nodelocs(6,:) = (loc_elem_nodelocs(2,:) + loc_elem_nodelocs(3,:)) ./ 2;
%    loc_elem_nodelocs(7,:) = (loc_elem_nodelocs(3,:) + loc_elem_nodelocs(1,:)) ./ 2;
%    loc_elem_nodelocs(8,:) = (loc_elem_nodelocs(1,:) + loc_elem_nodelocs(4,:)) ./ 2;
%    loc_elem_nodelocs(9,:) = (loc_elem_nodelocs(2,:) + loc_elem_nodelocs(4,:)) ./ 2;
%    loc_elem_nodelocs(10,:) = (loc_elem_nodelocs(3,:) + loc_elem_nodelocs(4,:)) ./ 2;


    % Modify local nodes (deformed configuration, using analytical deformation info)
    for m = 1:max(size(loc_elem_nodelocs))
      x = loc_elem_nodelocs(m,1);
      y = loc_elem_nodelocs(m,2);
      z = loc_elem_nodelocs(m,3);
      w = -1*x^2 + R^2;
      loc_elem_nodelocs(m,3) = z+w;
    end
  
    loc1_bary = zeros(1,4);
    loc1_bary(loc_elem_surf) = int_loc1;
    loc1 = tet_interp(loc_elem_nodelocs, loc1_bary);
    r1 = sqrt(loc1(1)^2 + loc1(2)^2);
    t1 = atan2(loc1(2), loc1(1));
    anal1_normal = circular_def_normal(loc1);
  
    loc2_bary = zeros(1,4);
    loc2_bary(loc_elem_surf) = int_loc2;
    loc2 = tet_interp(loc_elem_nodelocs, loc2_bary);
    r2 = sqrt(loc2(1)^2 + loc2(2)^2);
    t2 = atan2(loc2(2), loc2(1));
    anal2_normal = circular_def_normal(loc2);
  
    loc3_bary = zeros(1,4);
    loc3_bary(loc_elem_surf) = int_loc3;
    loc3 = tet_interp(loc_elem_nodelocs, loc3_bary);
    r3 = sqrt(loc3(1)^2 + loc3(2)^2);
    t3 = atan2(loc3(2), loc3(1));
    anal3_normal = circular_def_normal(loc3);
  
    % Determine numerically computed normal
    [num1_normal] = determine_normal(loc_elem_nodelocs, loc1, loc1_bary);
    [num2_normal] = determine_normal(loc_elem_nodelocs, loc2, loc2_bary);
    [num3_normal] = determine_normal(loc_elem_nodelocs, loc3, loc3_bary);
  
    diff1 = abs(norm(anal1_normal - num1_normal'));
    diff2 = abs(norm(anal2_normal - num2_normal'));
    diff3 = abs(norm(anal3_normal - num3_normal'));
    loc_diff = [diff1, diff2, diff3];
    diffs = [diffs, diff1, diff2, diff3];
    if (max(diffs) > max_diff)
      max_diff = max(diffs);
    end
  
    if ((max(loc_diff) > max_loc_diff) & (loc_elem_surf_num==2))
      max_loc_diff = max(loc_diff);
    end

    if (max_diff > 1.5)
      disp('stop'); pause;
    end
  
  end
  disp(strcat(['Max diff: ', num2str(max_diff)]));

  %diff_list = [diff_list, sqrt(sum(diffs))];
  diff_list = [diff_list, max_diff];
  loc_diff_list = [loc_diff_list, max_loc_diff];

end

c = diff_list(end) / elem_h_list(end);
yc = c.*elem_h_list;

c2 = diff_list(end) / elem_h_list(end)^2;
yc2 = c.*elem_h_list.^2;

%figure;
%loglog(elem_h_list, diff_list, '-ob'); hold on;
%loglog(elem_h_list, yc, '--r'); grid on;
%legend('Numerical', '$\mathcal{O} (h)$', 'interpreter', 'latex', 'fontsize', 18);
%xlabel('$h$', 'interpreter', 'latex', 'fontsize', 18);
%ylabel('$M$', 'interpreter', 'latex', 'fontsize', 18);
%title('Metric of Numerically Computed Unit Normals', 'interpreter', 'latex', 'fontsize', 18);

figure;
loglog(elem_h_list, diff_list, '-ob'); hold on;
loglog(elem_h_list, yc2, '--r'); grid on;
legend('Numerical', '$\mathcal{O} (h^2)$', 'interpreter', 'latex', 'fontsize', 18);
xlabel('$h$', 'interpreter', 'latex', 'fontsize', 18);
ylabel('$M$', 'interpreter', 'latex', 'fontsize', 18);
title('Metric of Numerically Computed Unit Normals', 'interpreter', 'latex', 'fontsize', 18);



function [normal_xyz] = circular_def_normal(location)

  loc_x = location(1);
  loc_y = location(2);
  loc_z = location(3);

  normal_xyz = [2*loc_x; 0; 1]; normal_xyz = normal_xyz ./ norm(normal_xyz);

end
