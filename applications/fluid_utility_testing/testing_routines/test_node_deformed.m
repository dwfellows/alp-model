
close all;
clear all;
clc;

fs = 18;

% Add mesh interfacing tools path
addpath('../../../fluid_function_utils');

% Pick random points and form their analytically derived surface normals
num_rand = 10;
analytical_surface_normals = zeros(3,num_rand);

R = 0.1;
k = 3.1962/R;
h = 0.0001;

% Pick random point
Abar = 0.001 / (besselj(0,0) - (besselj(0,k*R)/besseli(0,k*R))*besseli(0,0));
rand_points = zeros(2,num_rand);
for i = 1:num_rand
  % Define random point
  r = R*rand(1);
  t = 2*pi*rand(1);
  rand_points(:,i) = [r;t];

  % Define local surface normal analytically
  loc = [r*cos(t), r*sin(t), h];
  [anal_normal] = circular_def_normal(loc, Abar, R, k);
  analytical_surface_normals(:,i) = anal_normal;
end

% Work through random points on all meshes to evaluate convergence
for l = 1:5

  % Read in circular node data
  if (l==1)  % Coarse data
    elem_node = dlmread('../fem_data/circular_coarse/coarse.csv');
    node_elem = dlmread('../fem_data/circular_coarse/coarse_node_elem.csv');
    xdir_deform = dlmread('../fem_data/circular_coarse/coarse_xdir_deformation.txt');
    ydir_deform = dlmread('../fem_data/circular_coarse/coarse_ydir_deformation.txt');
    zdir_deform = dlmread('../fem_data/circular_coarse/coarse_zdir_deformation.txt');
    node_locs = xdir_deform(:,2:4);
    surface_node_data = dlmread('../fem_data/circular_coarse/coarse_surface.txt');
    coarse_num_normals = zeros(3,num_rand);
  elseif (l==2) %coarse-mid data
    elem_node = importdata('../fem_data/circular_cm/cm_elements.txt');
    elem_node = elem_node.data;
    node_elem = dlmread('../fem_data/circular_cm/cm_node_elem.csv');
    xdir_deform = dlmread('../fem_data/circular_cm/cm_xdir_deformation.txt');
    ydir_deform = dlmread('../fem_data/circular_cm/cm_ydir_deformation.txt');
    zdir_deform = dlmread('../fem_data/circular_cm/cm_zdir_deformation.txt');
    node_locs = xdir_deform(:,2:4);
    surface_node_data = dlmread('../fem_data/circular_cm/cm_surf_nodes.txt');
    cm_num_normals = zeros(3,num_rand);
  elseif (l==3)  % Mid data
    elem_node = dlmread('../fem_data/circular_mid/mid.csv');
    node_elem = dlmread('../fem_data/circular_mid/mid_node_elem.csv');
    xdir_deform = dlmread('../fem_data/circular_mid/mid_xdir_deformation.txt');
    ydir_deform = dlmread('../fem_data/circular_mid/mid_ydir_deformation.txt');
    zdir_deform = dlmread('../fem_data/circular_mid/mid_zdir_deformation.txt');
    node_locs = xdir_deform(:,2:4);
    surface_node_data = dlmread('../fem_data/circular_mid/mid_surface.txt');
    mid_num_normals = zeros(3,num_rand);
  elseif (l==4) % mid-fine data
    elem_node = importdata('../fem_data/circular_fm/fm_elements.txt');
    elem_node = elem_node.data;
    node_elem = dlmread('../fem_data/circular_fm/fm_node_elem.csv');
    xdir_deform = dlmread('../fem_data/circular_fm/fm_xdir_deformation.txt');
    ydir_deform = dlmread('../fem_data/circular_fm/fm_ydir_deformation.txt');
    zdir_deform = dlmread('../fem_data/circular_fm/fm_zdir_deformation.txt');
    node_locs = xdir_deform(:,2:4);
    surface_node_data = dlmread('../fem_data/circular_fm/fm_surface.txt');
    fm_num_normals = zeros(3,num_rand);
  else  % Fine data
    elem_node = dlmread('../fem_data/circular_fine/fine.csv');
    node_elem = dlmread('../fem_data/circular_fine/fine_node_elem.csv');
    xdir_deform = dlmread('../fem_data/circular_fine/fine_xdir_deformation.txt');
    ydir_deform = dlmread('../fem_data/circular_fine/fine_ydir_deformation.txt');
    zdir_deform = dlmread('../fem_data/circular_fine/fine_zdir_deformation.txt');
    node_locs = xdir_deform(:,2:4);
    surface_node_data = dlmread('../fem_data/circular_fine/fine_surface.txt');
    fine_num_normals = zeros(3,num_rand);
  end
  
  %% Create following single mode approximation
  % Pick random nodes, find their surface normals, compare to known analytical deformation of single mode approx
  N = max(size(surface_node_data)); 
  for i = 1:num_rand
  
    % Pick random point
    loc = rand_points(:,i);
    r = loc(1); t = loc(2);
    loc = [r*cos(t), r*sin(t), h];

    % Compute local surface normal numerically
    [local_struc_elem] = determine_elem(elem_node, node_elem, surface_node_data, node_locs, loc);
    local_struc_elem_nodes = elem_node(local_struc_elem, :);
    local_struc_elem_nodelocs = node_locs(local_struc_elem_nodes, :);
  
    % Compute local barycentric coordinates of test point (original configuration)
    loc_bary = cart2bary(local_struc_elem_nodelocs, loc);
  
    figure(1); clf;
    for m = 1:10
      scatter3(local_struc_elem_nodelocs(m,1), local_struc_elem_nodelocs(m,2), local_struc_elem_nodelocs(m,3), 'b'); hold on;
      text(local_struc_elem_nodelocs(m,1), local_struc_elem_nodelocs(m,2), local_struc_elem_nodelocs(m,3), num2str(m));
    end
    scatter3(loc(1), loc(2), loc(3), 'r');
    xlabel('$x$', 'interpreter', 'latex', 'fontsize', fs);
    ylabel('$y$', 'interpreter', 'latex', 'fontsize', fs);
    zlabel('$z$', 'interpreter', 'latex', 'fontsize', fs);
  
    % Modify local nodes (deformed configuration, using analytical deformation info)
    for m = 1:max(size(local_struc_elem_nodelocs))
      x = local_struc_elem_nodelocs(m,1);
      y = local_struc_elem_nodelocs(m,2);
      z = local_struc_elem_nodelocs(m,3);
      r = sqrt(x^2 + y^2);
      t = atan2(y,x);
      w = Abar*(besselj(0,k*r) - (besselj(0,k*R)/besseli(0,k*R))*besseli(0,k*r));
      local_struc_elem_nodelocs(m,3) = z+w;
    end
  
    % Interpolate test point deformation
    loc_def_z = tet_interp(local_struc_elem_nodelocs(:,3), loc_bary);
    %loc = [loc(1), loc(2), loc_def_z];
    loc = [tet_interp(local_struc_elem_nodelocs(:,1), loc_bary), tet_interp(local_struc_elem_nodelocs(:,2), loc_bary), tet_interp(local_struc_elem_nodelocs(:,3), loc_bary)];

    % Determine numerically computed normal
    [num_normal] = determine_normal(local_struc_elem_nodelocs, loc, loc_bary);
    if (l==1)
      coarse_num_normals(:,i) = num_normal;
    elseif (l==2)
      cm_num_normals(:,i) = num_normal;
    elseif (l==3)
      mid_num_normals(:,i) = num_normal;
    elseif (l==4)
      fm_num_normals(:,i) = num_normal;
    else
      fine_num_normals(:,i) = num_normal;
    end
  
    figure(2); clf;
    for m = 1:10
      scatter3(local_struc_elem_nodelocs(m,1), local_struc_elem_nodelocs(m,2), local_struc_elem_nodelocs(m,3), 'b'); hold on;
      text(local_struc_elem_nodelocs(m,1), local_struc_elem_nodelocs(m,2), local_struc_elem_nodelocs(m,3), num2str(m));
    end
    scatter3(loc(1), loc(2), loc(3), 'r');
    xlabel('$x$', 'interpreter', 'latex', 'fontsize', fs);
    ylabel('$y$', 'interpreter', 'latex', 'fontsize', fs);
    zlabel('$z$', 'interpreter', 'latex', 'fontsize', fs);
    quiver3(loc(1), loc(2), loc(3), anal_normal(1)/10000, anal_normal(2)/10000, anal_normal(3)/10000, 'r');
    quiver3(loc(1), loc(2), loc(3), num_normal(1)/10000, num_normal(2)/10000, num_normal(3)/10000, 'g');
  
  end
end

function [normal_xyz] = circular_def_normal(location, def_strength, radius, k)

  loc_x = location(1);
  loc_y = location(2);
  loc_z = location(3);
  Abar = def_strength;

  r = sqrt(loc_x^2 + loc_y^2);
  t = atan2(loc_y, loc_x);

  w = def_strength*(besselj(0,k*r) - (besselj(0,k*radius) / besseli(0,k*radius))*besseli(0,k*r));
  dwdr = def_strength*(-k*besselj(1,k*r) - k*(besselj(0,k*radius)/besseli(0,k*radius))*besseli(1,k*r));

  normal_xyz = [-dwdr*cos(t); -dwdr*sin(t); 1];
  normal_xyz = normal_xyz ./ norm(normal_xyz);

end
