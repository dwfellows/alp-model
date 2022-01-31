
close all;
clear all;
clc;

% This routine is used on an arbitrary geometry with defined shape function to
% test surface normal and surface normal derivative computations

fs = 18;

% Add mesh interfacing tools path
addpath('../../../fluid_function_utils');

% Pick random points to test at
num_rand = 1;
analytical_surface_normals = zeros(3, num_rand);
analytical_surface_normals_deriv = zeros(3, num_rand);

% Define attributes of circular plate
R = 0.1;
k = 3.1962/R;
h = 0.0001;

% Pick random points
Abar = 0.01 / (besselj(0,0) - (besselj(0,k*R)/besseli(0,k*R))*besseli(0,0));
rand_points = zeros(2,num_rand);
for i = 1:num_rand
  % Define random point
  r = R*rand(1);
  t = 2*pi*rand(1);
  r = 0.085820935324396899201992994221655;
  t = 2.7130624029126120255739351705415;
  rand_points(:,i) = [r;t];

  % Define local surface normal analytically
  loc = [r*cos(t), r*sin(t), h];
  [anal_normal, anal_normal_deriv] = circular_def_normal(loc, Abar, R, k);
  analytical_surface_normals(:,i) = anal_normal;
  analytical_surface_normals_deriv(:,i) = anal_normal_deriv;
end

% Work through random points on all meshes to evaluate convergence
for l = 1:5

  % Read in circular node data
  if (l==1) % Coarse data
    elem_node = dlmread('../fem_data/circular_coarse/coarse.csv');
    node_elem = dlmread('../fem_data/circular_coarse/coarse_node_elem.csv');
    xdir_deform = dlmread('../fem_data/circular_coarse/coarse_xdir_deformation.txt');
    ydir_deform = dlmread('../fem_data/circular_coarse/coarse_ydir_deformation.txt');
    zdir_deform = dlmread('../fem_data/circular_coarse/coarse_zdir_deformation.txt');
    node_locs = xdir_deform(:,2:4);
    surface_node_data = dlmread('../fem_data/circular_coarse/coarse_surface.txt');
    coarse_num_normals = zeros(3,num_rand);
    coarse_num_normals_deriv = zeros(3,num_rand);
  elseif (l==2) % cm data
    elem_node = importdata('../fem_data/circular_cm/cm_elements.txt');
    elem_node = elem_node.data;
    node_elem = dlmread('../fem_data/circular_cm/cm_node_elem.csv');
    xdir_deform = dlmread('../fem_data/circular_cm/cm_xdir_deformation.txt');
    ydir_deform = dlmread('../fem_data/circular_cm/cm_ydir_deformation.txt');
    zdir_deform = dlmread('../fem_data/circular_cm/cm_zdir_deformation.txt');
    node_locs = xdir_deform(:,2:4);
    surface_node_data = dlmread('../fem_data/circular_cm/cm_surf_nodes.txt');
    cm_num_normals = zeros(3,num_rand);
    cm_num_normals_deriv = zeros(3,num_rand);
  elseif (l==3)  % Mid data
    elem_node = dlmread('../fem_data/circular_mid/mid.csv');
    node_elem = dlmread('../fem_data/circular_mid/mid_node_elem.csv');
    xdir_deform = dlmread('../fem_data/circular_mid/mid_xdir_deformation.txt');
    ydir_deform = dlmread('../fem_data/circular_mid/mid_ydir_deformation.txt');
    zdir_deform = dlmread('../fem_data/circular_mid/mid_zdir_deformation.txt');
    node_locs = xdir_deform(:,2:4);
    surface_node_data = dlmread('../fem_data/circular_mid/mid_surface.txt');
    mid_num_normals = zeros(3,num_rand);
    mid_num_normals_deriv = zeros(3,num_rand);
  elseif (l==4) % fm data
    elem_node = importdata('../fem_data/circular_fm/fm_elements.txt');
    elem_node = elem_node.data;
    node_elem = dlmread('../fem_data/circular_fm/fm_node_elem.csv');
    xdir_deform = dlmread('../fem_data/circular_fm/fm_xdir_deformation.txt');
    ydir_deform = dlmread('../fem_data/circular_fm/fm_ydir_deformation.txt');
    zdir_deform = dlmread('../fem_data/circular_fm/fm_zdir_deformation.txt');
    node_locs = xdir_deform(:,2:4);
    surface_node_data = dlmread('../fem_data/circular_fm/fm_surface.txt');
    fm_num_normals = zeros(3,num_rand);
    fm_num_normals_deriv = zeros(3,num_rand);
  else  % Fine data
    elem_node = dlmread('../fem_data/circular_fine/fine.csv');
    node_elem = dlmread('../fem_data/circular_fine/fine_node_elem.csv');
    xdir_deform = dlmread('../fem_data/circular_fine/fine_xdir_deformation.txt');
    ydir_deform = dlmread('../fem_data/circular_fine/fine_ydir_deformation.txt');
    zdir_deform = dlmread('../fem_data/circular_fine/fine_zdir_deformation.txt');
    node_locs = xdir_deform(:,2:4);
    surface_node_data = dlmread('../fem_data/circular_fine/fine_surface.txt');
    fine_num_normals = zeros(3,num_rand);
    fine_num_normals_deriv = zeros(3,num_rand);
  end

  % Perform computational comparison
  N = max(size(surface_node_data));
  for i = 1:num_rand

    % Pick each random point
    loc = rand_points(:,i);
    r = loc(1); t = loc(2);
    z = (R^2 - r^2) + Abar*(besselj(0,k*r) - (besselj(0,k*R)/besseli(0,k*R))*besseli(0,k*r));
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
    for m = 1:max(size(local_struc_elem_nodelocs));
      x = local_struc_elem_nodelocs(m,1);
      y = local_struc_elem_nodelocs(m,2);
      z = local_struc_elem_nodelocs(m,3);
      r = sqrt(x^2 + y^2);
      t = atan2(y,x);
      local_struc_elem_nodelocs(m,3) = (R^2 - r^2) + Abar*(besselj(0,k*r) - (besselj(0,k*R)/besseli(0,k*R))*besseli(0,k*r));

    end
    
    % Interpolate test point deformation
    loc = [tet_interp(local_struc_elem_nodelocs(:,1), loc_bary), tet_interp(local_struc_elem_nodelocs(:,2), loc_bary), tet_interp(local_struc_elem_nodelocs(:,3), loc_bary)];

    % Construct deformation function at nodelocs
    def_shape_function = zeros(10,3);
    for m = 1:10;
      x = local_struc_elem_nodelocs(m,1);
      y = local_struc_elem_nodelocs(m,2);
      r = sqrt(x^2 + y^2);
      def_shape_function(m,:) = [0, 0, besselj(0,k*r) - (besselj(0,k*R)/besseli(0,k*R))*besseli(0,k*r)];
    end

    % Determine numerically computed normal and derivative
    [num_normal] = determine_normal(local_struc_elem_nodelocs, loc, loc_bary);
    [num_normal_deriv] = determine_normal_deriv(local_struc_elem_nodelocs, def_shape_function, loc, loc_bary);
    if (l==1)
      coarse_num_normals(:,i) = num_normal;
      coarse_num_normals_deriv(:,i) = num_normal_deriv;
    elseif (l==2)
      cm_num_normals(:,i) = num_normal;
      cm_num_normals_deriv(:,i) = num_normal_deriv;
    elseif (l==3)
      mid_num_normals(:,i) = num_normal;
      mid_num_normals_deriv(:,i) = num_normal_deriv;
    elseif (l==4)
      fm_num_normals(:,i) = num_normal;
      fm_num_normals_deriv(:,i) = num_normal_deriv;
    else
      fine_num_normals(:,i) = num_normal;
      fine_num_normals_deriv(:,i) = num_normal_deriv;
    end

  end
end

coarse_norm = vecnorm(coarse_num_normals_deriv - analytical_surface_normals_deriv) ./ vecnorm(analytical_surface_normals_deriv);
cm_norm = vecnorm(cm_num_normals_deriv - analytical_surface_normals_deriv) ./ vecnorm(analytical_surface_normals_deriv);
mid_norm = vecnorm(mid_num_normals_deriv - analytical_surface_normals_deriv) ./ vecnorm(analytical_surface_normals_deriv);
fm_norm = vecnorm(fm_num_normals_deriv - analytical_surface_normals_deriv) ./ vecnorm(analytical_surface_normals_deriv);
fine_norm = vecnorm(fine_num_normals_deriv - analytical_surface_normals_deriv) ./ vecnorm(analytical_surface_normals_deriv);
disp(strcat(['Coarse norm: ', num2str(coarse_norm)]));
disp(strcat(['CM norm: ', num2str(cm_norm)]));
disp(strcat(['Mid norm: ', num2str(mid_norm)]));
disp(strcat(['FM norm: ', num2str(fm_norm)]));
disp(strcat(['Fine norm: ', num2str(fine_norm)]));

function [normal_xyz, normal_deriv_xyz] = circular_def_normal(location, q, radius, k);

  loc_x = location(1);
  loc_y = location(2);
  loc_z = location(3);
  R = radius;

  r = sqrt(loc_x^2 + loc_y^2);
  t = atan2(loc_y, loc_x);

  m = (R^2 - r^2) + q*(besselj(0,k*r) - (besselj(0,k*R)/besseli(0,k*R))*besseli(0,k*r));
  dmdr = -2*r + q*(-k*besselj(1,k*r) - k*(besselj(0,k*R)/besseli(0,k*R))*besseli(1,k*r));
  d2mdrdq = -k*besselj(1,k*r) - k*(besselj(0,k*R)/besseli(0,k*R))*besseli(1,k*r);

  normal_xyz = [(-dmdr*cos(t))/(sqrt((dmdr^2)+1)); (-dmdr*sin(t))/(sqrt((dmdr^2)+1)); 1/sqrt((dmdr^2)+1)];

  normal_deriv_xyz = [(-d2mdrdq*cos(t)*(sqrt((dmdr^2)+1)-(dmdr^2)/sqrt((dmdr^2)+1)))/((dmdr^2)+1); ...
                      (-d2mdrdq*sin(t)*(sqrt((dmdr^2)+1)-(dmdr^2)/sqrt((dmdr^2)+1)))/((dmdr^2)+1); ...
                      (-(dmdr*d2mdrdq)/sqrt((dmdr^2)+1))/((dmdr^2)+1)];



end
