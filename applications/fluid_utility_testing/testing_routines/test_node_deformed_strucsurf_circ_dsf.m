
close all;
clear all;
clc;

fs = 24;

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

% Define maximum deformation
Abar = 0.001 / (besselj(0,0) - (besselj(0,k*R)/besseli(0,k*R))*besseli(0,0));

% Define integration point locations
int_loc1 = [2/3, 1/6, 1/6];
int_loc2 = [1/6, 2/3, 1/6];
int_loc3 = [1/6, 1/6, 2/3];

% Work through random points on all meshes to evaluate convergence
elem_h_list = [];
diff_list = []; loc_diff_list = [];
min_diff_list = [];
for l = [6 7 1:5]

  % Read in circular node data
  if (l==1)  % Coarse data
    elem_node = dlmread('../fem_data/circular_coarse/coarse.csv');
    node_elem = dlmread('../fem_data/circular_coarse/coarse_node_elem.csv');
    xdir_deform = dlmread('../fem_data/circular_coarse/coarse_xdir_deformation.txt');
    ydir_deform = dlmread('../fem_data/circular_coarse/coarse_ydir_deformation.txt');
    zdir_deform = dlmread('../fem_data/circular_coarse/coarse_zdir_deformation.txt');
    node_locs = xdir_deform(:,2:4);
    surface_node_data = dlmread('../fem_data/circular_coarse/coarse_surface.txt');
    struct_faces = dlmread('../fem_data/circular_coarse/struc_data/coarse_faces.txt');
    elem_a = area / max(size(struct_faces));
    elem_h_list = [elem_h_list, (elem_a*4/sqrt(3))^(1/2)];
  elseif (l==2) %coarse-mid data
    elem_node = importdata('../fem_data/circular_cm/cm_elements.txt');
    elem_node = elem_node.data;
    node_elem = dlmread('../fem_data/circular_cm/cm_node_elem.csv');
    xdir_deform = dlmread('../fem_data/circular_cm/cm_xdir_deformation.txt');
    ydir_deform = dlmread('../fem_data/circular_cm/cm_ydir_deformation.txt');
    zdir_deform = dlmread('../fem_data/circular_cm/cm_zdir_deformation.txt');
    node_locs = xdir_deform(:,2:4);
    surface_node_data = dlmread('../fem_data/circular_cm/cm_surf_nodes.txt');
    struct_faces = dlmread('../fem_data/circular_cm/struc_data/cm_faces.txt');
    elem_a = area / max(size(struct_faces));
    elem_h_list = [elem_h_list, (elem_a*4/sqrt(3))^(1/2)];
  elseif (l==3)  % Mid data
    elem_node = dlmread('../fem_data/circular_mid/mid.csv');
    node_elem = dlmread('../fem_data/circular_mid/mid_node_elem.csv');
    xdir_deform = dlmread('../fem_data/circular_mid/mid_xdir_deformation.txt');
    ydir_deform = dlmread('../fem_data/circular_mid/mid_ydir_deformation.txt');
    zdir_deform = dlmread('../fem_data/circular_mid/mid_zdir_deformation.txt');
    node_locs = xdir_deform(:,2:4);
    surface_node_data = dlmread('../fem_data/circular_mid/mid_surface.txt');
    struct_faces = dlmread('../fem_data/circular_mid/struc_data/mid_faces.txt');
    elem_a = area / max(size(struct_faces));
    elem_h_list = [elem_h_list, (elem_a*4/sqrt(3))^(1/2)];
  elseif (l==4) % mid-fine data
    elem_node = importdata('../fem_data/circular_fm/fm_elements.txt');
    elem_node = elem_node.data;
    node_elem = dlmread('../fem_data/circular_fm/fm_node_elem.csv');
    xdir_deform = dlmread('../fem_data/circular_fm/fm_xdir_deformation.txt');
    ydir_deform = dlmread('../fem_data/circular_fm/fm_ydir_deformation.txt');
    zdir_deform = dlmread('../fem_data/circular_fm/fm_zdir_deformation.txt');
    node_locs = xdir_deform(:,2:4);
    surface_node_data = dlmread('../fem_data/circular_fm/fm_surface.txt');
    struct_faces = dlmread('../fem_data/circular_fm/struc_data/fm_faces.txt');
    elem_a = area / max(size(struct_faces));
    elem_h_list = [elem_h_list, (elem_a*4/sqrt(3))^(1/2)];
  elseif (l==5)  % Fine data
    elem_node = dlmread('../fem_data/circular_fine/fine.csv');
    node_elem = dlmread('../fem_data/circular_fine/fine_node_elem.csv');
    xdir_deform = dlmread('../fem_data/circular_fine/fine_xdir_deformation.txt');
    ydir_deform = dlmread('../fem_data/circular_fine/fine_ydir_deformation.txt');
    zdir_deform = dlmread('../fem_data/circular_fine/fine_zdir_deformation.txt');
    node_locs = xdir_deform(:,2:4);
    surface_node_data = dlmread('../fem_data/circular_fine/fine_surface.txt');
    struct_faces = dlmread('../fem_data/circular_fine/struc_data/fine_faces.txt');
    elem_a = area / max(size(struct_faces));
    elem_h_list = [elem_h_list, (elem_a*4/sqrt(3))^(1/2)];
  elseif (l==6)
    elem_node = importdata('../fem_data/new_circ_data/circ_c1/circ_c1_elems.txt'); elem_node = elem_node.data;
    node_elem = dlmread('../fem_data/new_circ_data/circ_c1/circ_c1_node_elem.txt');
    xdir_deform = dlmread('../fem_data/new_circ_data/circ_c1/xdir_deform.txt');
    ydir_deform = dlmread('../fem_data/new_circ_data/circ_c1/ydir_deform.txt');
    zdir_deform = dlmread('../fem_data/new_circ_data/circ_c1/zdir_deform.txt');
    node_locs = dlmread('../fem_data/new_circ_data/circ_c1/circ_c1_nodes.txt'); node_locs = node_locs(:,2:4);
    surface_node_data = dlmread('../fem_data/new_circ_data/circ_c1/circ_c1_face_nodes.txt');
    struct_faces = dlmread('../fem_data/new_circ_data/circ_c1/circ_c1_faces.txt');
    elem_a = area / max(size(struct_faces));
    elem_h_list = [elem_h_list, (elem_a*4/sqrt(3))^(1/2)];
  elseif (l==7)
    elem_node = importdata('../fem_data/new_circ_data/circ_c2/circ_c2_elems.txt'); elem_node = elem_node.data;
    node_elem = dlmread('../fem_data/new_circ_data/circ_c2/circ_c2_node_elem.txt');
    xdir_deform = dlmread('../fem_data/new_circ_data/circ_c2/xdir_deform.txt');
    ydir_deform = dlmread('../fem_data/new_circ_data/circ_c2/ydir_deform.txt');
    zdir_deform = dlmread('../fem_data/new_circ_data/circ_c2/zdir_deform.txt');
    node_locs = dlmread('../fem_data/new_circ_data/circ_c2/circ_c2_nodes.txt'); node_locs = node_locs(:,2:4);
    surface_node_data = dlmread('../fem_data/new_circ_data/circ_c2/circ_c2_face_nodes.txt');
    struct_faces = dlmread('../fem_data/new_circ_data/circ_c2/circ_c2_faces.txt');
    elem_a = area / max(size(struct_faces));
    elem_h_list = [elem_h_list, (elem_a*4/sqrt(3))^(1/2)];
  end

  max_diff = 0; max_loc_diff = 0; min_diff = 2;
  counter = [0 0 0 0];
  diffs = [];
  location_diffs = [];

  def_shape_function = [xdir_deform(:,5), ydir_deform(:,5), zdir_deform(:,5)];
  max_def = 1000 * max( vecnorm(def_shape_function') ); 
  def_shape_function = def_shape_function ./ max_def;

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
    loc_dsf = def_shape_function(loc_elem_nodes,:);

    loc_elem_surf = [1 2 3 4]; loc_elem_surf(loc_elem_surf_num) = [];

    counter(loc_elem_surf_num) = counter(loc_elem_surf_num) + 1;
    
    % Modify local nodes to put mid-nodes on actual mid-line
%    loc_elem_nodelocs(5,:) = (loc_elem_nodelocs(1,:) + loc_elem_nodelocs(2,:)) ./ 2;
%    loc_elem_nodelocs(6,:) = (loc_elem_nodelocs(2,:) + loc_elem_nodelocs(3,:)) ./ 2;
%    loc_elem_nodelocs(7,:) = (loc_elem_nodelocs(3,:) + loc_elem_nodelocs(1,:)) ./ 2;
%    loc_elem_nodelocs(8,:) = (loc_elem_nodelocs(1,:) + loc_elem_nodelocs(4,:)) ./ 2;
%    loc_elem_nodelocs(9,:) = (loc_elem_nodelocs(2,:) + loc_elem_nodelocs(4,:)) ./ 2;
%    loc_elem_nodelocs(10,:) = (loc_elem_nodelocs(3,:) + loc_elem_nodelocs(4,:)) ./ 2;


    % Modify local nodes (deformed configuration, using deformation info)
    loc_elem_nodelocs = loc_elem_nodelocs + loc_dsf;
  
    loc1_bary = zeros(1,4);
    loc1_bary(loc_elem_surf) = int_loc1;
    loc1 = tet_interp(loc_elem_nodelocs, loc1_bary);
    r1 = sqrt(loc1(1)^2 + loc1(2)^2);
    t1 = atan2(loc1(2), loc1(1));
    anal1_normal = circular_def_normal(loc1, Abar, R, k, growth);
  
    loc2_bary = zeros(1,4);
    loc2_bary(loc_elem_surf) = int_loc2;
    loc2 = tet_interp(loc_elem_nodelocs, loc2_bary);
    r2 = sqrt(loc2(1)^2 + loc2(2)^2);
    t2 = atan2(loc2(2), loc2(1));
    anal2_normal = circular_def_normal(loc2, Abar, R, k, growth);
  
    loc3_bary = zeros(1,4);
    loc3_bary(loc_elem_surf) = int_loc3;
    loc3 = tet_interp(loc_elem_nodelocs, loc3_bary);
    r3 = sqrt(loc3(1)^2 + loc3(2)^2);
    t3 = atan2(loc3(2), loc3(1));
    anal3_normal = circular_def_normal(loc3, Abar, R, k, growth);
  
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
 
    if (min(diffs) < min_diff)
      min_diff = min(diffs);
    end
 
    if ((max(loc_diff) > max_loc_diff) & (loc_elem_surf_num==2))
      max_loc_diff = max(loc_diff);
    end

  end

  disp(strcat(['Max diff: ', num2str(max_diff)]));

  %diff_list = [diff_list, sqrt(sum(diffs))];
  diff_list = [diff_list, max_diff];
  loc_diff_list = [loc_diff_list, max_loc_diff];
  min_diff_list = [min_diff_list, min_diff];

end

c = diff_list(end) / elem_h_list(end);
yc = c.*elem_h_list;

%figure;
%loglog(elem_h_list, diff_list, '-ob'); hold on;
%loglog(elem_h_list, yc, '--r'); grid on;
%legend('Numerical', '$\mathcal{O} (h)$', 'interpreter', 'latex', 'location', 'northwest', 'fontsize', 18);
%xlabel('$h$', 'interpreter', 'latex', 'fontsize', 18);
%ylabel('$\max_{\Omega} \big( | \hat{n}_a - \hat{n}_n | \big)$', 'interpreter', 'latex', 'fontsize', 18);
%title('Metric of Numerically Computed Unit Normals using FEM Data', 'interpreter', 'latex', 'fontsize', 18);


c2 = diff_list(end) / elem_h_list(end)^2;
yc2 = c2.*elem_h_list.^2;

figure;
lg1 = loglog(elem_h_list, diff_list, '-ob'); hold on;
lg1.LineWidth = 2;
lg2 = loglog(elem_h_list, yc2, '--r'); grid on;
lg2.LineWidth = 2;
legend('Numerical', '$\mathcal{O} (h^2)$', 'interpreter', 'latex', 'location', 'northwest', 'fontsize', 24);
a = get(gca, 'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',18);
xlabel('Avg. Triangle Side Length [m]', 'interpreter', 'latex', 'fontsize', fs);
ylabel('$\max_{\Omega} \big( | \hat{n}_a - \hat{n}_n | \big)$', 'interpreter', 'latex', 'fontsize', fs);
title('$L_{\infty}$ Error of $\hat{n}$, Circular Plate, First Axisymmetric Mode (FEM Data)', 'interpreter', 'latex', 'fontsize', fs);



function [normal_xyz] = circular_def_normal(location, def_strength, radius, k, growth)

  loc_x = location(1);
  loc_y = location(2);
  loc_z = location(3);
  Abar = def_strength;

  r = sqrt(loc_x^2 + loc_y^2);
  t = atan2(loc_y, loc_x); t = wrapTo2Pi(t);

  w = def_strength*(besselj(0,k*r) - (besselj(0,k*radius) / besseli(0,k*radius))*besseli(0,k*r));
  dwdr = def_strength*(-k*besselj(1,k*r) - k*(besselj(0,k*radius)/besseli(0,k*radius))*besseli(1,k*r));

%  w = growth*(radius^2 - r^2);
%  dwdr = -2*growth*r;

  normal_xyz = [-dwdr*cos(t); -dwdr*sin(t); 1];
  normal_xyz = normal_xyz ./ norm(normal_xyz);

  %normal_xyz = [2*loc_x; 0; 1]; normal_xyz = normal_xyz ./ norm(normal_xyz);

end
