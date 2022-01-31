
close all;
clear all;
clc;

% This routine is used on an arbitrary geometry with defined shape function to
% test surface normal and surface normal derivative computations

fs = 24;

% Add mesh interfacing tools path
addpath('../../../fluid_function_utils');

% Define attributes of circular plate
R = 0.1;
besselfuns = @(l) besselj(0,l)*besseli(1,l) + besseli(0,l)*besselj(1,l);
lambda = fsolve(besselfuns, 3.1962);
k = lambda/R;
h = 0.0001;

area = pi*R^2;

% Pick random points
l1 = 1e-5;
Abar = l1 / (besselj(0,0) - (besselj(0,k*R)/besseli(0,k*R))*besseli(0,0));
l2 = (1/l1);

% Define integration point locations
int_loc1 = [2/3, 1/6, 1/6];
int_loc2 = [1/6, 2/3, 1/6];
int_loc3 = [1/6, 1/6, 2/3];

elem_h_list = []; max_diffs = [];
% Work through random points on all meshes to evaluate convergence
for l = [6 7 1:5]

  % Read in circular node data
  if (l==1) % Coarse data
    elem_node = dlmread('../fem_data/circular_coarse/coarse.csv');
    node_elem = dlmread('../fem_data/circular_coarse/coarse_node_elem.csv');
    xdir_deform = dlmread('../fem_data/circular_coarse/coarse_xdir_deformation.txt');
    ydir_deform = dlmread('../fem_data/circular_coarse/coarse_ydir_deformation.txt');
    zdir_deform = dlmread('../fem_data/circular_coarse/coarse_zdir_deformation.txt');
    node_locs = xdir_deform(:,2:4);
    surface_node_data = dlmread('../fem_data/circular_coarse/coarse_surface.txt');
    struct_faces = dlmread('../fem_data/circular_coarse/struc_data/coarse_faces.txt');
  elseif (l==2) % cm data
    elem_node = importdata('../fem_data/circular_cm/cm_elements.txt');
    elem_node = elem_node.data;
    node_elem = dlmread('../fem_data/circular_cm/cm_node_elem.csv');
    xdir_deform = dlmread('../fem_data/circular_cm/cm_xdir_deformation.txt');
    ydir_deform = dlmread('../fem_data/circular_cm/cm_ydir_deformation.txt');
    zdir_deform = dlmread('../fem_data/circular_cm/cm_zdir_deformation.txt');
    node_locs = xdir_deform(:,2:4);
    surface_node_data = dlmread('../fem_data/circular_cm/cm_surf_nodes.txt');
    struct_faces = dlmread('../fem_data/circular_cm/struc_data/cm_faces.txt');
  elseif (l==3)  % Mid data
    elem_node = dlmread('../fem_data/circular_mid/mid.csv');
    node_elem = dlmread('../fem_data/circular_mid/mid_node_elem.csv');
    xdir_deform = dlmread('../fem_data/circular_mid/mid_xdir_deformation.txt');
    ydir_deform = dlmread('../fem_data/circular_mid/mid_ydir_deformation.txt');
    zdir_deform = dlmread('../fem_data/circular_mid/mid_zdir_deformation.txt');
    node_locs = xdir_deform(:,2:4);
    surface_node_data = dlmread('../fem_data/circular_mid/mid_surface.txt');
    struct_faces = dlmread('../fem_data/circular_mid/struc_data/mid_faces.txt');
  elseif (l==4) % fm data
    elem_node = importdata('../fem_data/circular_fm/fm_elements.txt');
    elem_node = elem_node.data;
    node_elem = dlmread('../fem_data/circular_fm/fm_node_elem.csv');
    xdir_deform = dlmread('../fem_data/circular_fm/fm_xdir_deformation.txt');
    ydir_deform = dlmread('../fem_data/circular_fm/fm_ydir_deformation.txt');
    zdir_deform = dlmread('../fem_data/circular_fm/fm_zdir_deformation.txt');
    node_locs = xdir_deform(:,2:4);
    surface_node_data = dlmread('../fem_data/circular_fm/fm_surface.txt');
    struct_faces = dlmread('../fem_data/circular_fm/struc_data/fm_faces.txt');
  elseif (l==5)  % Fine data
    elem_node = dlmread('../fem_data/circular_fine/fine.csv');
    node_elem = dlmread('../fem_data/circular_fine/fine_node_elem.csv');
    xdir_deform = dlmread('../fem_data/circular_fine/fine_xdir_deformation.txt');
    ydir_deform = dlmread('../fem_data/circular_fine/fine_ydir_deformation.txt');
    zdir_deform = dlmread('../fem_data/circular_fine/fine_zdir_deformation.txt');
    node_locs = xdir_deform(:,2:4);
    surface_node_data = dlmread('../fem_data/circular_fine/fine_surface.txt');
    struct_faces = dlmread('../fem_data/circular_fine/struc_data/fine_faces.txt');
  elseif (l==6)
    elem_node = importdata('../fem_data/new_circ_data/circ_c1/circ_c1_elems.txt'); elem_node = elem_node.data;
    node_elem = dlmread('../fem_data/new_circ_data/circ_c1/circ_c1_node_elem.txt');
    xdir_deform = dlmread('../fem_data/new_circ_data/circ_c1/xdir_deform.txt');
    ydir_deform = dlmread('../fem_data/new_circ_data/circ_c1/ydir_deform.txt');
    zdir_deform = dlmread('../fem_data/new_circ_data/circ_c1/zdir_deform.txt');
    node_locs = dlmread('../fem_data/new_circ_data/circ_c1/circ_c1_nodes.txt'); node_locs = node_locs(:,2:4);
    surface_node_data = dlmread('../fem_data/new_circ_data/circ_c1/circ_c1_face_nodes.txt');
    struct_faces = dlmread('../fem_data/new_circ_data/circ_c1/circ_c1_faces.txt');
  elseif (l==7)
    elem_node = importdata('../fem_data/new_circ_data/circ_c2/circ_c2_elems.txt'); elem_node = elem_node.data;
    node_elem = dlmread('../fem_data/new_circ_data/circ_c2/circ_c2_node_elem.txt');
    xdir_deform = dlmread('../fem_data/new_circ_data/circ_c2/xdir_deform.txt');
    ydir_deform = dlmread('../fem_data/new_circ_data/circ_c2/ydir_deform.txt');
    zdir_deform = dlmread('../fem_data/new_circ_data/circ_c2/zdir_deform.txt');
    node_locs = dlmread('../fem_data/new_circ_data/circ_c2/circ_c2_nodes.txt'); node_locs = node_locs(:,2:4);
    surface_node_data = dlmread('../fem_data/new_circ_data/circ_c2/circ_c2_face_nodes.txt');
    struct_faces = dlmread('../fem_data/new_circ_data/circ_c2/circ_c2_faces.txt');
  end

  elem_a = area / max(size(struct_faces));
  elem_h_list = [elem_h_list, sqrt(elem_a*(4/sqrt(3)))];

  max_diff_norm = 0;
  max_diff_nd = 0;
  max_anal_norm = 0;

  % If using numerical deformation data
  def_shape_function = [xdir_deform(:,5), ydir_deform(:,5), zdir_deform(:,5)];
  def_norms = vecnorm(def_shape_function'); def_norms_max = l2 * max(def_norms);
  def_shape_function = def_shape_function ./ (def_norms_max);

  % Perform computational comparisons
  N = max(size(struct_faces));
  for i = 1:N

    % Construct integration points
    loc_elem = struct_faces(i,1);
    loc_elem_surf_nodes = struct_faces(i, 2:11);

    loc_elem_nodes = elem_node(loc_elem, :);
    loc_elem_nodelocs = node_locs(loc_elem_nodes, :);
    loc_elem_surf_num = find(loc_elem_surf_nodes(1:4) == 0);
    loc_elem_surf = [1 2 3 4]; loc_elem_surf(loc_elem_surf_num) = [];

    % Construct deformation function and pre-deformation at nodelocs analytically
    for m = 1:10;
      x = loc_elem_nodelocs(m,1);
      y = loc_elem_nodelocs(m,2);
      r = sqrt(x^2 + y^2);
      loc_def_shape_function(m,:) = [0, 0, Abar*(besselj(0,k*r) - (besselj(0,k*R)/besseli(0,k*R))*besseli(0,k*r))];
      z = loc_elem_nodelocs(m,3);
    end

    % Use FEM prescribed deformation function
    loc_def_shape_function = def_shape_function(loc_elem_nodes, :);

    r = sqrt(loc_elem_nodelocs(:,1).^2 + loc_elem_nodelocs(:,2).^2);
    t = atan2(loc_elem_nodelocs(:,2), loc_elem_nodelocs(:,1));

    loc1_bary = zeros(1,4);
    loc1_bary(loc_elem_surf) = int_loc1;
    loc1 = tet_interp(loc_elem_nodelocs, loc1_bary);
    r1 = sqrt(loc1(1)^2 + loc1(2)^2);
    t1 = atan2(loc1(2), loc1(1));
    [anal1_normal, anal1_normal_deriv] = circular_def_normal(loc1, Abar, R, k);

    loc2_bary = zeros(1,4);
    loc2_bary(loc_elem_surf) = int_loc2;
    loc2 = tet_interp(loc_elem_nodelocs, loc2_bary);
    r2 = sqrt(loc2(1)^2 + loc2(2)^2);
    t2 = atan2(loc2(2), loc2(1));
    [anal2_normal, anal2_normal_deriv] = circular_def_normal(loc2, Abar, R, k);

    loc3_bary = zeros(1,4);
    loc3_bary(loc_elem_surf) = int_loc3;
    loc3 = tet_interp(loc_elem_nodelocs, loc3_bary);
    r3 = sqrt(loc3(1)^2 + loc3(2)^2);
    t3 = atan2(loc3(2), loc3(1));
    [anal3_normal, anal3_normal_deriv] = circular_def_normal(loc3, Abar, R, k);

    % Determine numerically computed normal and derivative
    [num1_normal] = determine_normal(loc_elem_nodelocs, loc1, loc1_bary);
    [num2_normal] = determine_normal(loc_elem_nodelocs, loc2, loc2_bary);
    [num3_normal] = determine_normal(loc_elem_nodelocs, loc3, loc3_bary);

    [num1_normal_deriv] = determine_normal_deriv(loc_elem_nodelocs, loc_def_shape_function, loc1, loc1_bary);
    [num2_normal_deriv] = determine_normal_deriv(loc_elem_nodelocs, loc_def_shape_function, loc2, loc2_bary);
    [num3_normal_deriv] = determine_normal_deriv(loc_elem_nodelocs, loc_def_shape_function, loc3, loc3_bary);

    diff1n = abs(norm(num1_normal' - anal1_normal));
    diff2n = abs(norm(num2_normal' - anal2_normal));
    diff3n = abs(norm(num3_normal' - anal3_normal));

    diff1nd = abs(norm(num1_normal_deriv' - anal1_normal_deriv));
    diff2nd = abs(norm(num2_normal_deriv' - anal2_normal_deriv));
    diff3nd = abs(norm(num3_normal_deriv' - anal3_normal_deriv));

    anal1norm = norm(anal1_normal_deriv);
    anal2norm = norm(anal2_normal_deriv);
    anal3norm = norm(anal3_normal_deriv);
    analnorms = [anal1norm, anal2norm, anal3norm];

    diffsn = [diff1n, diff2n, diff3n];
    diffsnd = [diff1nd, diff2nd, diff3nd];

    if (max(diffsn) > max_diff_norm)
      max_diff_norm = max(diffsn);
    end

    if (max(diffsnd) > max_diff_nd)
      max_diff_nd = max(diffsnd);
      max_diff_elem = loc_elem;
    end

    if (max(analnorms) > max_anal_norm)
      max_anal_norm = max(analnorms);
    end

  end

  disp(strcat(['Max normal diff: ', num2str(max_diff_norm)]));
  disp(strcat(['Max normal deriv diff: ', num2str(max_diff_nd)]));

  % Correct norm to take, after discussion with Alex. Divide max error by max correct vector.
  max_diffs = [max_diffs, max_diff_nd / max_anal_norm];

end

c = max_diffs(end) / elem_h_list(end);
yc = c.*elem_h_list;

%figure; hold on; grid on;
%loglog(elem_h_list, max_diffs, '-ob');
%loglog(elem_h_list, yc, '--r');
%legend('Numerical', '$\mathcal{O}(h)$', 'interpreter', 'latex', 'location', 'northwest', 'fontsize', fs);
%xlabel('$h$', 'interpreter', 'latex', 'fontsize', fs);
%ylabel('$\Big( \max_{\Omega} || \frac{d\hat{n}}{dq_1} |_{n} - \frac{d\hat{n}}{dq_1} |_{a} ||  \Big) / \Big( \max_{\Omega} || \frac{d\hat{n}}{dq_1} |_a || \Big) $', 'interpreter', 'latex', 'fontsize', 18);
%title('Matrix of Errors in Derivative of Surface Unit Normal, FEM Data', 'interpreter', 'latex', 'fontsize', fs);

c2 = max_diffs(end) / elem_h_list(end)^2;
yc2 = c2.*elem_h_list.^2;

figure;
lg1 = loglog(elem_h_list, max_diffs, '-ob'); hold on;
lg1.LineWidth = 2;
lg2 = loglog(elem_h_list, yc2, '--r'); grid on;
lg2.LineWidth = 2;
a = get(gca, 'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',18);
legend('Numerical', '$\mathcal{O}(h^2)$', 'interpreter', 'latex', 'location', 'northwest', 'fontsize', 24);
xlabel('Avg. Triangle Side Length [m]', 'interpreter', 'latex', 'fontsize', fs);
ylabel('$\Big( \max_{\Omega} || \frac{d\hat{n}}{dq_1} |_{n} - \frac{d\hat{n}}{dq_1} |_{a} ||  \Big) / \Big( \max_{\Omega} || \frac{d\hat{n}}{dq_1} |_a || \Big) $', 'interpreter', 'latex', 'fontsize', fs);
title('$L_{\infty}$ Error of $\frac{d\hat{n}}{dq_1}$, Circular Plate, First Axisymmetric Mode (FEM Data)', 'interpreter', 'latex', 'fontsize', fs);


function [normal_xyz, normal_deriv_xyz] = circular_def_normal(location, Abar, radius, k);

  loc_x = location(1);
  loc_y = location(2);
  loc_z = location(3);
  R = radius;

  r = sqrt(loc_x^2 + loc_y^2);
  t = atan2(loc_y, loc_x);  t = wrapTo2Pi(t);

  m = 0;
  dmdr = 0;
  d2mdrdq = Abar*(-k*besselj(1,k*r) - k*(besselj(0,k*R)/besseli(0,k*R))*besseli(1,k*r));

  normal_xyz = [(-dmdr*cos(t))/(sqrt((dmdr^2)+1)); (-dmdr*sin(t))/(sqrt((dmdr^2)+1)); 1/sqrt((dmdr^2)+1)];

%   normal_deriv_xyz = [(-d2mdrdq*cos(t)*(sqrt((dmdr^2)+1)-(dmdr^2)/sqrt((dmdr^2)+1)))/((dmdr^2)+1); ...
%                       (-d2mdrdq*sin(t)*(sqrt((dmdr^2)+1)-(dmdr^2)/sqrt((dmdr^2)+1)))/((dmdr^2)+1); ...
%                       (-(dmdr*d2mdrdq)/sqrt((dmdr^2)+1))/((dmdr^2)+1)];

  normal_deriv_xyz = zeros(3,1);
  normal_deriv_xyz(1) = (1/(dmdr^2+1))*( sqrt(dmdr^2 + 1)*(-d2mdrdq*cos(t)) - (-dmdr*cos(t))*((dmdr*d2mdrdq)/sqrt(dmdr^2 + 1)));
  normal_deriv_xyz(2) = (1/(dmdr^2+1))*( sqrt(dmdr^2 + 1)*(-d2mdrdq*sin(t)) - (-dmdr*sin(t))*((dmdr*d2mdrdq)/sqrt(dmdr^2 + 1)));
  normal_deriv_xyz(3) = (-dmdr*d2mdrdq)/(((dmdr)^2+1)^(3/2));

end

function [returned_locs] = correct_locs(given_locs)

    returned_locs = zeros(10,3);
    returned_locs(1:4,:) = given_locs(1:4,:);

    returned_locs(5,:) = (given_locs(1,:) + given_locs(2,:)) ./ 2;
    returned_locs(6,:) = (given_locs(2,:) + given_locs(3,:)) ./ 2;
    returned_locs(7,:) = (given_locs(3,:) + given_locs(1,:)) ./ 2;
    returned_locs(8,:) = (given_locs(1,:) + given_locs(4,:)) ./ 2;
    returned_locs(9,:) = (given_locs(2,:) + given_locs(4,:)) ./ 2;
    returned_locs(10,:) = (given_locs(3,:) + given_locs(4,:)) ./ 2;

end
