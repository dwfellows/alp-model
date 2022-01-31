
% Perform quadrature about structural element faces

close all;
clear all;
clc;

addpath('../../../fluid_function_utils');

fs = 26;

%% Import fluid data -- all quantities, from surface
% 1 node num | 2 X [m] | 3 Y [m] | 4 Z [m] | 5 p_a [Pa] | 6 A [m^2] | 7 \rho [kg m^-3] | 8 F_x [N] | 9 F_y [N] 
% 10 F_z [N] | 11 a [m s^-1] | 12 Mach # | 13 Mach # Stn | 14 normal_x | 15 normal_y | 16 normal_z | 17 v_u stn [m s^1]
% 18 v_v stn [m s^-1] | 19 v_w stn [m s^-1] | 20 v_u [m s^-1] | 21 v_v [m s^-1] | 22 v_w [m s^-1] | 23 X [m]
% 24 Y [m] | 25 Z [m]
%fluid_data = dlmread('../turbine_blade_data/free_slip/rotating_inner_wall_fluid_data.csv');
%face_conn = dlmread('../turbine_blade_data/free_slip/rotating_inner_wall_face_conn.csv');
%fluid_rho = fluid_data(:,7);
%fluid_p = fluid_data(:,5);
%fluid_a = fluid_data(:,11);
%fluid_vu = fluid_data(:,20);
%fluid_vv = fluid_data(:,21);
%fluid_vw = fluid_data(:,22);
%fluid_node_locs = fluid_data(:,2:4);

% Construct fluid data using coarse nodes -- to understand impact of fluid mesh error
%fluid_node_surf_rect = dlmread('../circular_data/new_circ_data/rect_fluid/rect_fluid_data.csv');
%face_conn = dlmread('../circular_data/new_circ_data/rect_fluid/rect_face_conn.csv');
%fluid_node_locs = fluid_node_surf_rect(:,2:4);
%fluid_node_faces = dlmread('../circular_data/new_circ_data/rect_fluid/rect_node_faces.txt');

% Construct fluid data using fine nodes
fluid_node_surf_rect = dlmread('../fem_data/new_circ_data/rect_fluid_hires/rect_fluid_data.csv');
face_conn = dlmread('../fem_data/new_circ_data/rect_fluid_hires/rect_face_conn.csv');
fluid_node_locs = fluid_node_surf_rect(:,2:4);
fluid_node_faces = dlmread('../fem_data/new_circ_data/rect_fluid_hires/rect_node_faces.txt');


% Define square plate parameters
L = 0.2; h = 0.0002;
area = L^2;

% Construct deformation function
%l1 = 0.0001; l2 = 1/l1;
l1 = 0.00001; l2 = 1/l1;
%l1 = 1; l2 = 1/l1;

abar = l1 ./ ((cos(pi)-1)^2)
psi3_anal = @(x,y) abar.*( (cos(2*pi.*x./L) - 1).*(cos(2*pi.*y./L) - 1) );

% Define analytical values
%rhol = @(x,y) ((10.*x./L).^6).*((10.*y./L).^4) + ((10.*x./L).^3) + ((10.*y./L).^2) + 10;
%rhol = @(x,y) (-y.*(y-L)) .* x.^1 + 1; % Works
%rhol = @(x,y) (-1.*y.*(y-L)) .* x.^0 + 1; % Works
%rhol = @(x,y) (-1.*x.*(x-L)) .* y.^0; % Works
rhol = @(x,y) (-1.*x.*(x-L)) .* y.^1;  % Problematic -- likely just an artefact of the local mesh
%rhol = @(x,y) 1.*(x.^1).*(y.^0);  % Works
%rhol = @(x,y) 1.*(x.^0).*(y.^0);
%rhol = @(x,y) 1 + x + y;  % Works
%rhol = @(x,y) (x.^0).*y.^1 + 1;
%rhol = @(x,y) (y.^0).*(x.^1) + 1;
al = @(x,y) sqrt(1.4*287*300) .* ( x.^0 ) .* (y.^0);
vl = @(x,y) [1;1;0] .* (x.^0) .* (y.^0);

% Compute surace normal quantities analytically
dwdx = @(x,y) abar .* ( (cos(2*pi.*y./L)-1).*( (-2*pi./L).*sin(2*pi.*x./L) ) );
dwdy = @(x,y) abar .* ( (cos(2*pi.*x./L)-1).*( (-2*pi./L).*sin(2*pi.*y./L) ) );
n1_anal = @(x,y) 0.*(x.^0).*(y.^0);
n2_anal = @(x,y) 0.*(x.^0).*(y.^0);
n3_anal = @(x,y) 1*(x.^0).*(y.^0);
dn1dq_anal = @(x,y) -1.*dwdx(x,y);
dn2dq_anal = @(x,y) -1.*dwdy(x,y);
dn3dq_anal = @(x,y) 0.*(x.^0).*(y.^0);


% Define quadrature locations (integration points) -- utilize 3-pt rule as FEM data is second-order accurate
% Rule from Alex and is contained in FEM textbooks
int_loc1 = [2/3, 1/6, 1/6];
int_loc2 = [1/6, 2/3, 1/6];
int_loc3 = [1/6, 1/6, 2/3];

dF1dq_approx = [];
n_err = [];
dn_err = [];
dn_2err = [];
ip_loc = [];
elem_h_list = [];
face_num = [];
mean_loc_dndq_err = [];
max_dndq_err_loc = [];

%  deformation_errs = zeros(2,7);
%  for l = 1:7
%    if (l==1)
%      elem_node = importdata('../circular_data/rect_plate_thin/rect_c1/elements.txt'); elem_node = elem_node.data;
%      node_elem = dlmread('../circular_data/rect_plate_thin/rect_c1/node_elem.txt');
%      nodes = dlmread('../circular_data/rect_plate_thin/rect_c1/nodes.txt');
%      structural_node_locs = nodes(:,2:4);
%      xdir_deform = dlmread('../circular_data/rect_plate_thin/rect_c1/xdir_deform.txt');
%      ydir_deform = dlmread('../circular_data/rect_plate_thin/rect_c1/ydir_deform.txt');
%      zdir_deform = dlmread('../circular_data/rect_plate_thin/rect_c1/zdir_deform.txt');
%      surface_node_data = dlmread('../circular_data/rect_plate_thin/rect_c1/face_nodes.txt');
%      struct_faces = dlmread('../circular_data/rect_plate_thin/rect_c1/surf_faces.txt');
%    elseif (l==2)
%      elem_node = importdata('../circular_data/rect_plate_thin/rect_c2/elements.txt'); elem_node = elem_node.data;
%      node_elem = dlmread('../circular_data/rect_plate_thin/rect_c2/node_elem.txt');
%      nodes = dlmread('../circular_data/rect_plate_thin/rect_c2/nodes.txt');
%      structural_node_locs = nodes(:,2:4);
%      xdir_deform = dlmread('../circular_data/rect_plate_thin/rect_c2/xdir_deform.txt');
%      ydir_deform = dlmread('../circular_data/rect_plate_thin/rect_c2/ydir_deform.txt');
%      zdir_deform = dlmread('../circular_data/rect_plate_thin/rect_c2/zdir_deform.txt');
%      surface_node_data = dlmread('../circular_data/rect_plate_thin/rect_c2/face_nodes.txt');
%      struct_faces = dlmread('../circular_data/rect_plate_thin/rect_c2/surf_faces.txt');
%    elseif (l==3)
%      elem_node = importdata('../circular_data/rect_plate_thin/rect_c3/elements.txt'); elem_node = elem_node.data;
%      node_elem = dlmread('../circular_data/rect_plate_thin/rect_c3/node_elem.txt');
%      nodes = dlmread('../circular_data/rect_plate_thin/rect_c3/nodes.txt');
%      structural_node_locs = nodes(:,2:4);
%      xdir_deform = dlmread('../circular_data/rect_plate_thin/rect_c3/xdir_deform.txt');
%      ydir_deform = dlmread('../circular_data/rect_plate_thin/rect_c3/ydir_deform.txt');
%      zdir_deform = dlmread('../circular_data/rect_plate_thin/rect_c3/zdir_deform.txt');
%      surface_node_data = dlmread('../circular_data/rect_plate_thin/rect_c3/face_nodes.txt');
%      struct_faces = dlmread('../circular_data/rect_plate_thin/rect_c3/surf_faces.txt');
%    elseif (l==4)
%      elem_node = importdata('../circular_data/rect_plate_thin/rect_c4/elements.txt'); elem_node = elem_node.data;
%      node_elem = dlmread('../circular_data/rect_plate_thin/rect_c4/node_elem.txt');
%      nodes = dlmread('../circular_data/rect_plate_thin/rect_c4/nodes.txt');
%      structural_node_locs = nodes(:,2:4);
%      xdir_deform = dlmread('../circular_data/rect_plate_thin/rect_c4/xdir_deform.txt');
%      ydir_deform = dlmread('../circular_data/rect_plate_thin/rect_c4/ydir_deform.txt');
%      zdir_deform = dlmread('../circular_data/rect_plate_thin/rect_c4/zdir_deform.txt');
%      surface_node_data = dlmread('../circular_data/rect_plate_thin/rect_c4/face_nodes.txt');
%      struct_faces = dlmread('../circular_data/rect_plate_thin/rect_c4/surf_faces.txt');
%    elseif (l==5)
%      elem_node = importdata('../circular_data/rect_plate_thin/rect_c5/elements.txt'); elem_node = elem_node.data;
%      node_elem = dlmread('../circular_data/rect_plate_thin/rect_c5/node_elem.txt');
%      nodes = dlmread('../circular_data/rect_plate_thin/rect_c5/nodes.txt');
%      structural_node_locs = nodes(:,2:4);
%      xdir_deform = dlmread('../circular_data/rect_plate_thin/rect_c5/xdir_deform.txt');
%      ydir_deform = dlmread('../circular_data/rect_plate_thin/rect_c5/ydir_deform.txt');
%      zdir_deform = dlmread('../circular_data/rect_plate_thin/rect_c5/zdir_deform.txt');
%      surface_node_data = dlmread('../circular_data/rect_plate_thin/rect_c5/face_nodes.txt');
%      struct_faces = dlmread('../circular_data/rect_plate_thin/rect_c5/surf_faces.txt');
%    elseif (l==6)
%      elem_node = importdata('../circular_data/rect_plate_thin/rect_c6/elements.txt'); elem_node = elem_node.data;
%      node_elem = dlmread('../circular_data/rect_plate_thin/rect_c6/node_elem.txt');
%      nodes = dlmread('../circular_data/rect_plate_thin/rect_c6/nodes.txt');
%      structural_node_locs = nodes(:,2:4);
%      xdir_deform = dlmread('../circular_data/rect_plate_thin/rect_c6/xdir_deform.txt');
%      ydir_deform = dlmread('../circular_data/rect_plate_thin/rect_c6/ydir_deform.txt');
%      zdir_deform = dlmread('../circular_data/rect_plate_thin/rect_c6/zdir_deform.txt');
%      surface_node_data = dlmread('../circular_data/rect_plate_thin/rect_c6/face_nodes.txt');
%      struct_faces = dlmread('../circular_data/rect_plate_thin/rect_c6/surf_faces.txt');
%    elseif (l==7)
%      elem_node = importdata('../circular_data/rect_plate_thin/rect_c7/elements.txt'); elem_node = elem_node.data;
%      node_elem = dlmread('../circular_data/rect_plate_thin/rect_c7/node_elem.txt');
%      nodes = dlmread('../circular_data/rect_plate_thin/rect_c7/nodes.txt');
%      structural_node_locs = nodes(:,2:4);
%      xdir_deform = dlmread('../circular_data/rect_plate_thin/rect_c7/xdir_deform.txt');
%      ydir_deform = dlmread('../circular_data/rect_plate_thin/rect_c7/ydir_deform.txt');
%      zdir_deform = dlmread('../circular_data/rect_plate_thin/rect_c7/zdir_deform.txt');
%      surface_node_data = dlmread('../circular_data/rect_plate_thin/rect_c7/face_nodes.txt');
%      struct_faces = dlmread('../circular_data/rect_plate_thin/rect_c7/surf_faces.txt');
%    end
%  
%    % Compute characteristic size of element
%    N = max(size(struct_faces));
%    elem_a = area / N;
%    elem_h_list = [elem_h_list, sqrt(elem_a*4/sqrt(3)) / L];
%  
%    % Construct FEM data
%    def_shape_function = [xdir_deform(:,5), ydir_deform(:,5), zdir_deform(:,5)];
%    def_shape_function_norms = vecnorm(def_shape_function'); max_def = l2*max(def_shape_function_norms);
%    def_shape_function = def_shape_function ./ max_def;
%  
%    psi_err = zeros(max(size(surface_node_data)), 3);
%    for i = 1:max(size(surface_node_data))
%      node = surface_node_data(i,1);
%      xl = surface_node_data(i,2); yl = surface_node_data(i,3);
%      psi_err(i,:) = [xl, yl, abs(def_shape_function(node,3) - psi3_anal(xl,yl))];
%    end
%    figure;
%    scatter3(psi_err(:,1), psi_err(:,2), psi_err(:,3));
%    xlabel('$x$', 'interpreter', 'latex', 'fontsize', fs); 
%    ylabel('$y$', 'interpreter', 'latex', 'fontsize', fs); 
%    zlabel('$z$', 'interpreter', 'latex', 'fontsize', fs);
%    title(strcat(['Local Error, FEM Solution and Analytical Solution, $h/L = ', num2str(elem_h_list(end)), '$, Thin Square Plate']), 'interpreter', 'latex', 'fontsize', fs);
%    
%  
%    deformation_errs(:,l) = [max(psi_err(:,3)); max(def_shape_function(:,1))/max(def_shape_function(:,3))];
%  end
%  
%  c2 = deformation_errs(1,end) / elem_h_list(end)^2;
%  yc2 = c2.*elem_h_list.^2;
%  
%  figure;
%  lg1 = loglog(elem_h_list, deformation_errs(1,:), '-ob'); hold on;
%  lg1.LineWidth = 2;
%  lg2 = loglog(elem_h_list, yc2, '--r'); grid on;
%  lg2.LineWidth = 2;
%  legend('Numerical', '$\mathcal{O}(h^2)$', 'interpreter', 'latex', 'fontsize', 18);
%  xlabel('$h/L$', 'interpreter', 'latex', 'fontsize', fs); 
%  ylabel('$\max_{\Omega} \Big( \psi_{3} - w_{3} \Big)$', 'interpreter', 'latex', 'fontsize', fs);
%  title('$L_{\infty}$ Error of FEM Solution to Thin Plate Approx.', 'interpreter', 'latex', 'fontsize', fs);

elem_h_list = [];

%% Import structural data -- all quantities
for l = 1:7
  if (l==1)
    elem_node = importdata('../fem_data/rect_plate_thin/rect_c1/elements.txt'); elem_node = elem_node.data;
    node_elem = dlmread('../fem_data/rect_plate_thin/rect_c1/node_elem.txt');
    nodes = dlmread('../fem_data/rect_plate_thin/rect_c1/nodes.txt');
    structural_node_locs = nodes(:,2:4);
    xdir_deform = dlmread('../fem_data/rect_plate_thin/rect_c1/xdir_deform.txt');
    ydir_deform = dlmread('../fem_data/rect_plate_thin/rect_c1/ydir_deform.txt');
    zdir_deform = dlmread('../fem_data/rect_plate_thin/rect_c1/zdir_deform.txt');
    surface_node_data = dlmread('../fem_data/rect_plate_thin/rect_c1/face_nodes.txt');
    struct_faces = dlmread('../fem_data/rect_plate_thin/rect_c1/surf_faces.txt');
  elseif (l==2)
    elem_node = importdata('../fem_data/rect_plate_thin/rect_c2/elements.txt'); elem_node = elem_node.data;
    node_elem = dlmread('../fem_data/rect_plate_thin/rect_c2/node_elem.txt');
    nodes = dlmread('../fem_data/rect_plate_thin/rect_c2/nodes.txt');
    structural_node_locs = nodes(:,2:4);
    xdir_deform = dlmread('../fem_data/rect_plate_thin/rect_c2/xdir_deform.txt');
    ydir_deform = dlmread('../fem_data/rect_plate_thin/rect_c2/ydir_deform.txt');
    zdir_deform = dlmread('../fem_data/rect_plate_thin/rect_c2/zdir_deform.txt');
    surface_node_data = dlmread('../fem_data/rect_plate_thin/rect_c2/face_nodes.txt');
    struct_faces = dlmread('../fem_data/rect_plate_thin/rect_c2/surf_faces.txt');
  elseif (l==3)
    elem_node = importdata('../fem_data/rect_plate_thin/rect_c3/elements.txt'); elem_node = elem_node.data;
    node_elem = dlmread('../fem_data/rect_plate_thin/rect_c3/node_elem.txt');
    nodes = dlmread('../fem_data/rect_plate_thin/rect_c3/nodes.txt');
    structural_node_locs = nodes(:,2:4);
    xdir_deform = dlmread('../fem_data/rect_plate_thin/rect_c3/xdir_deform.txt');
    ydir_deform = dlmread('../fem_data/rect_plate_thin/rect_c3/ydir_deform.txt');
    zdir_deform = dlmread('../fem_data/rect_plate_thin/rect_c3/zdir_deform.txt');
    surface_node_data = dlmread('../fem_data/rect_plate_thin/rect_c3/face_nodes.txt');
    struct_faces = dlmread('../fem_data/rect_plate_thin/rect_c3/surf_faces.txt');
  elseif (l==4)
    elem_node = importdata('../fem_data/rect_plate_thin/rect_c4/elements.txt'); elem_node = elem_node.data;
    node_elem = dlmread('../fem_data/rect_plate_thin/rect_c4/node_elem.txt');
    nodes = dlmread('../fem_data/rect_plate_thin/rect_c4/nodes.txt');
    structural_node_locs = nodes(:,2:4);
    xdir_deform = dlmread('../fem_data/rect_plate_thin/rect_c4/xdir_deform.txt');
    ydir_deform = dlmread('../fem_data/rect_plate_thin/rect_c4/ydir_deform.txt');
    zdir_deform = dlmread('../fem_data/rect_plate_thin/rect_c4/zdir_deform.txt');
    surface_node_data = dlmread('../fem_data/rect_plate_thin/rect_c4/face_nodes.txt');
    struct_faces = dlmread('../fem_data/rect_plate_thin/rect_c4/surf_faces.txt');
  elseif (l==5)
    elem_node = importdata('../fem_data/rect_plate_thin/rect_c5/elements.txt'); elem_node = elem_node.data;
    node_elem = dlmread('../fem_data/rect_plate_thin/rect_c5/node_elem.txt');
    nodes = dlmread('../fem_data/rect_plate_thin/rect_c5/nodes.txt');
    structural_node_locs = nodes(:,2:4);
    xdir_deform = dlmread('../fem_data/rect_plate_thin/rect_c5/xdir_deform.txt');
    ydir_deform = dlmread('../fem_data/rect_plate_thin/rect_c5/ydir_deform.txt');
    zdir_deform = dlmread('../fem_data/rect_plate_thin/rect_c5/zdir_deform.txt');
    surface_node_data = dlmread('../fem_data/rect_plate_thin/rect_c5/face_nodes.txt');
    struct_faces = dlmread('../fem_data/rect_plate_thin/rect_c5/surf_faces.txt');
  elseif (l==6)
    elem_node = importdata('../fem_data/rect_plate_thin/rect_c6/elements.txt'); elem_node = elem_node.data;
    node_elem = dlmread('../fem_data/rect_plate_thin/rect_c6/node_elem.txt');
    nodes = dlmread('../fem_data/rect_plate_thin/rect_c6/nodes.txt');
    structural_node_locs = nodes(:,2:4);
    xdir_deform = dlmread('../fem_data/rect_plate_thin/rect_c6/xdir_deform.txt');
    ydir_deform = dlmread('../fem_data/rect_plate_thin/rect_c6/ydir_deform.txt');
    zdir_deform = dlmread('../fem_data/rect_plate_thin/rect_c6/zdir_deform.txt');
    surface_node_data = dlmread('../fem_data/rect_plate_thin/rect_c6/face_nodes.txt');
    struct_faces = dlmread('../fem_data/rect_plate_thin/rect_c6/surf_faces.txt');
  elseif (l==7)
    elem_node = importdata('../fem_data/rect_plate_thin/rect_c7/elements.txt'); elem_node = elem_node.data;
    node_elem = dlmread('../fem_data/rect_plate_thin/rect_c7/node_elem.txt');
    nodes = dlmread('../fem_data/rect_plate_thin/rect_c7/nodes.txt');
    structural_node_locs = nodes(:,2:4);
    xdir_deform = dlmread('../fem_data/rect_plate_thin/rect_c7/xdir_deform.txt');
    ydir_deform = dlmread('../fem_data/rect_plate_thin/rect_c7/ydir_deform.txt');
    zdir_deform = dlmread('../fem_data/rect_plate_thin/rect_c7/zdir_deform.txt');
    surface_node_data = dlmread('../fem_data/rect_plate_thin/rect_c7/face_nodes.txt');
    struct_faces = dlmread('../fem_data/rect_plate_thin/rect_c7/surf_faces.txt');
  end

  % Construct FEM data
  def_shape_function = [xdir_deform(:,5), ydir_deform(:,5), zdir_deform(:,5)];
  def_shape_function_norms = vecnorm(def_shape_function'); max_def = l2*max(def_shape_function_norms);
  def_shape_function = def_shape_function ./ max_def;

  for i = 1:max(size(node_elem))
    xl = xdir_deform(i,2); yl = xdir_deform(i,3);
    def_shape_function(i,:) = [0 0 psi3_anal(xl,yl)];
  end

%  psi_err = zeros(max(size(surface_node_data)), 3);
%  for i = 1:max(size(surface_node_data))
%    node = surface_node_data(i,1);
%    xl = surface_node_data(i,2); yl = surface_node_data(i,3);
%    psi_err(i,:) = [xl, yl, abs(def_shape_function(node,3) - psi3_anal(xl,yl))];
%  end
%  figure;
%  scatter3(psi_err(:,1), psi_err(:,2), psi_err(:,3));
%  xlabel('x'); ylabel('y'); zlabel('z');
%  disp(max(psi_err(:,3))); pause;   

  % Rotate node locations 90-deg
%  R = -1.*[0 -1 0; 1 0 0; 0 0 0];
%  structural_node_locs = R*(structural_node_locs');
%  structural_node_locs = structural_node_locs';
%  structural_node_locs(:,2) = structural_node_locs(:,2) + L;

  %% Construct flow functions on fluid surface
  N = max(size(fluid_node_locs));
  fluid_rho = zeros(N,1); fluid_a = zeros(N,1);
  fluid_vu = zeros(N,1); fluid_vv = zeros(N,1); fluid_vw = zeros(N,1);
  for j = 1:max(size(fluid_node_locs))
    xl = fluid_node_locs(j,1); yl = fluid_node_locs(j,2);
    fluid_rho(j) = rhol(xl,yl);
    fluid_a(j) = al(xl,yl);
    temp_v = vl(xl,yl);
    fluid_vu(j) = temp_v(1); fluid_vv(j) = temp_v(2);
    fluid_vw(j) = temp_v(3);
  end


  % Compute characteristic size of element
  N = max(size(struct_faces));
  elem_a = area / N;
  elem_h_list = [elem_h_list, sqrt(elem_a*4/sqrt(3))];
  face_num = [face_num, N];

  % Keep local n error and dndq error
  loc_n_err = [];
  loc_dndq_err = [];
  loc_psi_err = [];
  loc_dndq = [];
  dndq_err_locs = [];

  %% Compute dF1dq on local mesh
  disp(strcat(['Mesh Case Number: ', num2str(l)]));
  disp(strcat(['Number of faces: ', num2str(N)]));
  dF1dq = 0;
  for i = 1:N
    disp(strcat(['Face ', num2str(i)]));
  
    % Grab face element number and nodes on surface
    loc_elem = struct_faces(i,1);
    loc_elem_surf_nodes = struct_faces(i,2:11);

    % Get all element nodes and locations
    loc_elem_nodes = elem_node(loc_elem, :);
    loc_elem_nodelocs = structural_node_locs(loc_elem_nodes,:);
    loc_elem_surf_num = find(loc_elem_surf_nodes(1:4) == 0);
    loc_elem_surf = [1 2 3 4]; loc_elem_surf(loc_elem_surf_num) = [];

    %%% Go over each int loc and form integrand
    %%% Int loc 1
    % Get location of integration point in Cartesian (global) and barycentric (local) coords
    location1_bary = zeros(1,4); 
    location1_bary(loc_elem_surf) = int_loc1;
    location1 = tet_interp(loc_elem_nodelocs, location1_bary);
  
    % Interpolate deformation, normal and normal derivative at integration point
    psi_nodes = def_shape_function(loc_elem_nodes,:);
    psi1 = tet_interp(psi_nodes, location1_bary);
    n1 = determine_normal(loc_elem_nodelocs, location1, location1_bary, loc_elem_surf_num);
    dn1dq = determine_normal_deriv(loc_elem_nodelocs, psi_nodes, location1, location1_bary, loc_elem_surf_num);
    [dxdLi1, dxdLj1, ~] = def_derivs(loc_elem_nodelocs, psi_nodes, location1, location1_bary, loc_elem_surf_num);
    A1 = cross(dxdLi1, dxdLj1); A1 = norm(A1);
  
%    % Determine fluid face containing integration point and face coordinates of structural point in the face
%    [loc1_fluid_face, loc1_struct_bary] = determine_fluid_face(face_conn, fluid_node_surf_rect, fluid_node_faces, location1);
%    loc1_fluid_face_nodes = face_conn(loc1_fluid_face, :) + 1;
%  
%    % Construct local fluid face nodes for interpolation of rho, a, velocity
%    loc1_rho_nodes = fluid_rho(loc1_fluid_face_nodes);
%    loc1_a_nodes = fluid_a(loc1_fluid_face_nodes);
%    loc1_vu_nodes = fluid_vu(loc1_fluid_face_nodes);
%    loc1_vv_nodes = fluid_vv(loc1_fluid_face_nodes);
%    loc1_vw_nodes = fluid_vw(loc1_fluid_face_nodes);
%  
%    % Perform triangular interpolation
%    rho1 = tri_interp(loc1_struct_bary, loc1_rho_nodes);
%    a1 = tri_interp(loc1_struct_bary, loc1_a_nodes);
%    vu1 = tri_interp(loc1_struct_bary, loc1_vu_nodes);
%    vv1 = tri_interp(loc1_struct_bary, loc1_vv_nodes);
%    vw1 = tri_interp(loc1_struct_bary, loc1_vw_nodes);
%    v1 = [vu1; vv1; vw1];
  
    xl = location1(1); yl = location1(2);
    rho1 = rhol(xl,yl);
    a1 = al(xl,yl);
    vu1 = 1; vv1 = 1; vw1 = 0; v1 = [vu1; vv1; vw1];
    loc_dndq_err = [loc_dndq_err, norm(dn1dq'-norm_deriv_dq(xl,yl,L,abar))]; dndq_err_locs = [dndq_err_locs, i];
    loc_dndq = [loc_dndq, norm(norm_deriv_dq(xl,yl,L,abar))];
    ip_loc = [ip_loc, [xl;yl]];
    dn1dq = norm_deriv_dq(xl,yl,L,abar); dn1dq = dn1dq';
    loc_psi_err = [loc_psi_err, norm(psi1' - [0;0;psi3_anal(xl,yl)])];

    %%% Int loc 2
    % Get location of integration point in Cartesian (global) and barycentric (local) coords
    location2_bary = zeros(1,4);
    location2_bary(loc_elem_surf) = int_loc2;
    location2 = tet_interp(loc_elem_nodelocs, location2_bary);
  
    % Interpolate deformation, normal, and normal derivative at integration point
    psi2 = tet_interp(psi_nodes, location2_bary);
    n2 = determine_normal(loc_elem_nodelocs, location2, location2_bary, loc_elem_surf_num);
    dn2dq = determine_normal_deriv(loc_elem_nodelocs, psi_nodes, location2, location2_bary, loc_elem_surf_num);
    [dxdLi2, dxdLj2, ~] = def_derivs(loc_elem_nodelocs, psi_nodes, location2, location2_bary, loc_elem_surf_num);
    A2 = cross(dxdLi2, dxdLj2); A2 = norm(A2);
  
    % Determine fluid face containing integration point and face coordinates of structural point in the face
%    [loc2_fluid_face, loc2_struct_bary] = determine_fluid_face(face_conn, fluid_node_surf_rect, fluid_node_faces, location2);
%    loc2_fluid_face_nodes = face_conn(loc2_fluid_face, :) + 1;
%  
%    % Construct local fluid face nodes for interpolation of rho, a, velocity
%    loc2_rho_nodes = fluid_rho(loc2_fluid_face_nodes);
%    loc2_a_nodes = fluid_a(loc2_fluid_face_nodes);
%    loc2_vu_nodes = fluid_vu(loc2_fluid_face_nodes);
%    loc2_vv_nodes = fluid_vv(loc2_fluid_face_nodes);
%    loc2_vw_nodes = fluid_vw(loc2_fluid_face_nodes);
%  
%    % Perform triangular interpolation
%    rho2 = tri_interp(loc2_struct_bary, loc2_rho_nodes);
%    a2 = tri_interp(loc2_struct_bary, loc2_a_nodes);
%    vu2 = tri_interp(loc2_struct_bary, loc2_vu_nodes);
%    vv2 = tri_interp(loc2_struct_bary, loc2_vv_nodes);
%    vw2 = tri_interp(loc2_struct_bary, loc2_vw_nodes);
%    v2 = [vu2; vv2; vw2];

    xl = location2(1); yl = location2(2);
    rho2 = rhol(xl,yl);
    a2 = al(xl,yl);
    vu2 = 1; vv2 = 1; vw2 = 0; v2 = [vu2; vv2; vw2];
    loc_dndq_err = [loc_dndq_err, norm(dn2dq' - norm_deriv_dq(xl,yl,L,abar))]; dndq_err_locs = [dndq_err_locs, i];
    loc_dndq = [loc_dndq, norm(norm_deriv_dq(xl,yl,L,abar))];
    ip_loc = [ip_loc, [xl;yl]];
    dn2dq = norm_deriv_dq(xl,yl,L,abar); dn2dq = dn2dq';
    loc_psi_err = [loc_psi_err, norm(psi2' - [0;0;psi3_anal(xl,yl)])];

    %%% Int loc 3
    % Get location of integration point in Cartesian (global) and barycentric (local) coords
    location3_bary = zeros(1,4);
    location3_bary(loc_elem_surf) = int_loc3;
    location3 = tet_interp(loc_elem_nodelocs, location3_bary);

    % Interpolate deformation, normal, and normal derivative at integration point
    psi3 = tet_interp(psi_nodes, location3_bary);
    n3 = determine_normal(loc_elem_nodelocs, location3, location3_bary, loc_elem_surf_num);
    dn3dq = determine_normal_deriv(loc_elem_nodelocs, psi_nodes, location3, location3_bary, loc_elem_surf_num);
    [dxdLi3, dxdLj3, ~] = def_derivs(loc_elem_nodelocs, psi_nodes, location3, location3_bary, loc_elem_surf_num);
    A3 = cross(dxdLi3, dxdLj3); A3 = norm(A3);
  
%    % Determine fluid face containing integration point and face coordinates of structural point in the face
%    [loc3_fluid_face, loc3_struct_bary] = determine_fluid_face(face_conn, fluid_node_surf_rect, fluid_node_faces, location3);
%    loc3_fluid_face_nodes = face_conn(loc3_fluid_face, :) + 1;
%  
%    % Construct local fluid face nodes for interpolation of rho, a, velocity
%    loc3_rho_nodes = fluid_rho(loc3_fluid_face_nodes);
%    loc3_a_nodes = fluid_a(loc3_fluid_face_nodes);
%    loc3_vu_nodes = fluid_vu(loc3_fluid_face_nodes);
%    loc3_vv_nodes = fluid_vv(loc3_fluid_face_nodes);
%    loc3_vw_nodes = fluid_vw(loc3_fluid_face_nodes);
%  
%    % Perform triangular interpolation
%    rho3 = tri_interp(loc3_struct_bary, loc3_rho_nodes);
%    a3 = tri_interp(loc3_struct_bary, loc3_a_nodes);
%    vu3 = tri_interp(loc3_struct_bary, loc3_vu_nodes);
%    vv3 = tri_interp(loc3_struct_bary, loc3_vv_nodes);
%    vw3 = tri_interp(loc3_struct_bary, loc3_vw_nodes);
%    v3 = [vu3; vv3; vw3];

    xl = location3(1); yl = location3(2);
    rho3 = rhol(xl,yl);
    a3 = al(xl,yl);
    vu3 = 1; vv3 = 1; vw3 = 0; v3 = [vu3; vv3; vw3];
    loc_dndq_err = [loc_dndq_err, norm(dn3dq' - norm_deriv_dq(xl,yl,L,abar))]; dndq_err_locs = [dndq_err_locs, i];
    loc_dndq = [loc_dndq, norm(norm_deriv_dq(xl,yl,L,abar))];
    ip_loc = [ip_loc, [xl;yl]];
    dn3dq = norm_deriv_dq(xl,yl,L,abar); dn3dq = dn3dq';
    loc_psi_err = [loc_psi_err, norm(psi3' - [0;0;psi3_anal(xl,yl)])];

    if ((norm(n2 - [0 0 1])) > 1e-8)
      disp('hey'); pause;
    end
  
    % Compute integrands at integration points
    dF1dq_int1 = -rho1*a1*(dot(v1,dn1dq')*dot(n1',psi1') + dot(v1,n1')*dot(dn1dq',psi1'))*A1;
    dF1dq_int2 = -rho2*a2*(dot(v2,dn2dq')*dot(n2',psi2') + dot(v2,n2')*dot(dn2dq',psi2'))*A2;
    dF1dq_int3 = -rho3*a3*(dot(v3,dn3dq')*dot(n3',psi3') + dot(v3,n3')*dot(dn3dq',psi3'))*A3;
  
    % Perform quadrature
    dF1dq_loc = quadrature(dF1dq_int1, dF1dq_int2, dF1dq_int3);
    dF1dq = dF1dq + dF1dq_loc;

    loc_n_err = [loc_n_err, norm([0;0;1]-n1'), norm([0;0;1]-n2'), norm([0;0;1]-n3')];

  end

  dF1dq_approx = [dF1dq_approx, dF1dq];
  n_err = [n_err, max(loc_n_err)];
  dn_err = [dn_err, max(loc_dndq_err) / max(loc_dndq)];
  dn_2err = [dn_2err, sqrt(sum((loc_dndq_err ./ max(loc_dndq)).^2))];

%  figure;
%  scatter3(ip_loc(1,:), ip_loc(2,:), loc_dndq_err ./ max(loc_dndq));
%  xlabel('x'); ylabel('y'); zlabel('z');
%  ip_loc = [];

%  figure;
%  scatter3(ip_loc(1,:), ip_loc(2,:), loc_psi_err);
%  xlabel('x'); ylabel('y'); zlabel('z');
%  title(strcat(['Local $\vec{\psi}$ Error, Mesh ', num2str(l)]));

  mean_loc_dndq_err = [mean_loc_dndq_err, mean(loc_dndq_err)];
  [max_dndq_err, max_dndq_err_ind] = max(loc_dndq_err);
  max_dndq_err_loc = [max_dndq_err_loc, dndq_err_locs(max_dndq_err_ind)];

end

% Compute solution analytically
fun = @(x,y) rhol(x,y) .* al(x,y) .* ((dn1dq_anal(x,y) + dn2dq_anal(x,y)).*(n3_anal(x,y).*psi3_anal(x,y)) + (n1_anal(x,y) + n2_anal(x,y)).*(dn3dq_anal(x,y).*psi3_anal(x,y)));
dF1dq_anal = -1 * integral2(fun, 0, L, 0, L);

% Compute difference metrics

% Non-Dimensionalize Triangle Side Lengths
elem_h_list = elem_h_list ./ L;

if (abs(dF1dq_anal) > 1e-10)  % Compute absolute and relative errors

  % Absolute
  diffs = abs(dF1dq_approx - dF1dq_anal);

  c2 = diffs(end) / elem_h_list(end)^2;
  yc2 = c2.*elem_h_list.^2;

  figure;
  lg1 = loglog(elem_h_list, diffs, '-ob'); hold on;
  lg1.LineWidth = 2;
  lg2 = loglog(elem_h_list, yc2, '--r'); grid on;
  lg2.LineWidth = 2;
  legend('Numerical', '$\mathcal{O}(h^2)$', 'interpreter', 'latex', 'location', 'northwest', 'fontsize', 24);
  a = get(gca,'YTickLabel');
  set(gca,'YTickLabel',a,'fontsize',18);
  grid on;
  xlabel('$h / L$', 'interpreter', 'latex', 'fontsize', fs);
  ylabel('Absolute Error: $|| I_n - I_a ||$', 'interpreter', 'latex', 'fontsize', fs);
  title('Numerical Approximation of Imposed $\frac{dF_1}{dq_1} (q_1^e = 0)$', 'interpreter', 'latex', 'fontsize', fs);

  % Relative
  diffs = abs(dF1dq_approx - dF1dq_anal) ./ abs(dF1dq_anal);

  c2 = diffs(end) / elem_h_list(end)^2;
  yc2 = c2.*elem_h_list.^2;

  figure;
  lg1 = loglog(elem_h_list, diffs, '-ob'); hold on;
  lg1.LineWidth = 2;
  lg2 = loglog(elem_h_list, yc2, '--r'); grid on;
  lg2.LineWidth = 2;
  legend('Numerical', '$\mathcal{O}(h^2)$', 'interpreter', 'latex', 'location', 'northwest', 'fontsize', 24);
  a = get(gca,'YTickLabel');
  set(gca,'YTickLabel',a,'fontsize',18);
  grid on;
  xlabel('$h / L$', 'interpreter', 'latex', 'fontsize', fs);
  ylabel('Relative Error: $|| I_n - I_a || / || I_a ||$', 'interpreter', 'latex', 'fontsize', fs);
  title('Numerical Approximation of Imposed $\frac{dF_1}{dq_1} (q_1^e = 0)$', 'interpreter', 'latex', 'fontsize', fs);

else  % Compute only absolute errors

  diffs = abs(dF1dq_approx - dF1dq_anal);

  c2 = diffs(end) / elem_h_list(end)^2;
  yc2 = c2.*elem_h_list.^2;

  figure;
  lg1 = loglog(elem_h_list, diffs, '-ob'); hold on;
  lg1.LineWidth = 2;
  lg2 = loglog(elem_h_list, yc2, '--r'); grid on;
  lg2.LineWidth = 2;
  legend('Numerical', '$\mathcal{O}(h^2)$', 'interpreter', 'latex', 'location', 'northwest', 'fontsize', 24);
  a = get(gca,'YTickLabel');
  set(gca,'YTickLabel',a,'fontsize',18);
  grid on;
  xlabel('$h / L$', 'interpreter', 'latex', 'fontsize', fs);
  ylabel('Absolute Error: $|| I_n - I_a ||$', 'interpreter', 'latex', 'fontsize', fs);
  title('Numerical Approximation of Imposed $\frac{dF_1}{dq_1} (q_1^e = 0)$', 'interpreter', 'latex', 'fontsize', fs);

end

% Plot average local error on log-log plot
diffs = abs(dF1dq_approx - dF1dq_anal) ./ face_num;
c2 = diffs(end) / elem_h_list(end)^2;
yc2 = c2.*elem_h_list.^2;

figure;
lg1 = loglog(elem_h_list, diffs, '-ob'); hold on;
lg1.LineWidth = 2;
lg2 = loglog(elem_h_list, yc2, '--r'); grid on;
lg2.LineWidth = 2;
legend('Numerical', '$\mathcal{O}(h^2)$', 'interpreter', 'latex', 'location', 'northwest', 'fontsize', 24);
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',18);
grid on;
xlabel('$h / L$', 'interpreter', 'latex', 'fontsize', fs);
ylabel('Local Error: $|| I_n - I_a || / N$', 'interpreter', 'latex', 'fontsize', fs);
title('Local Error from Numerical Approximation of Imposed $\frac{dF_1}{dq_1} (q_1^e = 0)$', 'interpreter', 'latex', 'fontsize', fs);


% Convergence to final solution
diffs = abs(dF1dq_approx(1:end-1) - dF1dq_approx(end)) ./ abs(dF1dq_approx(end));
c2 = diffs(end) / elem_h_list(end-1)^2;
yc2 = c2.*elem_h_list(1:end-1).^2;

figure;
lg1 = loglog(elem_h_list(1:end-1), diffs, '-ob'); hold on;
lg1.LineWidth = 2;
lg2 = loglog(elem_h_list(1:end-1), yc2, '--r'); grid on;
lg2.LineWidth = 2;
legend('Numerical', '$\mathcal{O}(h^2)$', 'interpreter', 'latex', 'location', 'northwest', 'fontsize', 24);
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',18);
grid on;
xlabel('$h / L$', 'interpreter', 'latex', 'fontsize', fs);
ylabel('Max Error: $|| \max_{\Omega} \frac{d\hat{n}}{dq_1} |_{n} - \frac{d\hat{n}}{dq_1} |_{n,f} || / || \max_{\Omega} \frac{d\hat{n}}{dq_1} |_{n,f} ||$', 'interpreter', 'latex', 'fontsize', fs);
title('$L_{\infty}$ Error of Self-Convergence from Numerical Approximation of $\frac{d\hat{n}}{dq_1}$', 'interpreter', 'latex', 'fontsize', fs);


% Plot error metric of dn/dq
c2 = dn_err(end) / elem_h_list(end)^2;
yc2 = c2.*elem_h_list.^2;

figure;
lg1 = loglog(elem_h_list, dn_err, '-ob'); hold on;
lg1.LineWidth = 2;
lg2 = loglog(elem_h_list, yc2, '--r'); grid on;
lg2.LineWidth = 2;
legend('Numerical', '$\mathcal{O}(h^2)$', 'interpreter', 'latex', 'location', 'northwest', 'fontsize', 24);
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',18);
grid on;
xlabel('$h / L$', 'interpreter', 'latex', 'fontsize', fs);
ylabel('Max Error: $|| \max_{\Omega} \frac{d\hat{n}}{dq_1} |_{n} - \frac{d\hat{n}}{dq_1} |_{a} || / || \max_{\Omega} \frac{d\hat{n}}{dq_1} |_{a} ||$', 'interpreter', 'latex', 'fontsize', fs);
title('$L_{\infty}$ Error from Numerical Approximation of $\frac{d\hat{n}}{dq_1}$', 'interpreter', 'latex', 'fontsize', fs);


% Plot the L2 error metric of dn/dq
c2 = dn_2err(end) / elem_h_list(end)^2;
yc2 = c2.*elem_h_list.^2;

figure;
lg1 = loglog(elem_h_list, dn_2err, '-ob'); hold on;
lg1.LineWidth = 2;
lg2 = loglog(elem_h_list, yc2, '--r'); grid on;
lg2.LineWidth = 2;
legend('Numerical', '$\mathcal{O}(h^2)$', 'interpreter', 'latex', 'location', 'northwest', 'fontsize', 24);
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',18);
grid on;
xlabel('$h / L$', 'interpreter', 'latex', 'fontsize', fs);
ylabel('$L_2$ Error: $\Bigg( \sum || \frac{d\hat{n}}{dq_1} |_{n} - \frac{d\hat{n}}{dq_1} |_{a} || / \max_{\Omega} \frac{d\hat{n}}{dq_1} |_{a} || \Bigg)^{1/2}$', 'interpreter', 'latex', 'fontsize', fs);
title('$L_2$ Error from Numerical Approximation of $\frac{d\hat{n}}{dq_1}$', 'interpreter', 'latex', 'fontsize', fs);


%  % Compare self convergence
%  if (abs(dF1dq_approx(end)) > 1e-10)
%    diffs = abs(dF1dq_approx(1:end-1) - dF1dq_approx(end)) / abs(dF1dq_approx(end));
%  else
%    diffs = abs(dF1dq_approx(1:end-1) - dF1dq_approx(end));
%  end
%  
%  c = diffs(end) / elem_h_list(end-1);
%  yc = c.*elem_h_list(1:end-1);
%  
%  c2 = diffs(end) / (elem_h_list(end-1)^2);
%  yc2 = c2.*elem_h_list(1:end-1).^2;
%  
%  figure;
%  lg1 = loglog(elem_h_list(1:end-1), diffs, '-ob'); hold on;
%  lg1.LineWidth = 2;
%  lg2 = loglog(elem_h_list(1:end-1), yc2, '--r'); grid on;
%  lg2.LineWidth = 2;
%  legend('Numerical', '$\mathcal{O}(h^2)$', 'interpreter', 'latex', 'location', 'northwest', 'fontsize', 24);
%  a = get(gca,'YTickLabel');
%  set(gca,'YTickLabel',a,'fontsize',18);
%  grid on;
%  xlabel('$Avg. Triangle Side Length [m]$', 'interpreter', 'latex', 'fontsize', fs);
%  ylabel('Error: $|| I_n - I_{n,f} || / || I_{n,f} ||$', 'interpreter', 'latex', 'fontsize', fs);
%  title('FEM Self Convergence of Approximation of Imposed $\frac{dF_1}{dq_1} (q_1^e = 0)$', 'interpreter', 'latex', 'fontsize', fs);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Triangular interpolation function for location
function [loc] = tri_loc(face_node_locs, interp_loc)
  % Linear triangle location from area coordinates
  % face_node_locs: Cartesian locations of each face node
  % interp_loc: location of point in triangle, area coordinates
  % loc: location of point in triangle, Cartesian coordinates

  loc = face_node_locs(1,:).*interp_loc(1) + face_node_locs(2,:).*interp_loc(2) + face_node_locs(3,:).*interp_loc(3);

end

% Triangular interpolation function for data
function [quant] = tri_interp(interp_loc, data);

  quant = data(1)*interp_loc(1) + data(2)*interp_loc(2) + data(3)*interp_loc(3);

end

% Third-order quadrature function
function [int] = quadrature(val1, val2, val3)

  % Weights obtained from Zienkiewics and modified to ensure that transformed surface integral provides A = 1/2.
  w1 = 1/6;
  w2 = 1/6;
  w3 = 1/6;

  int = w1*val1 + w2*val2 + w3*val3;

end


function [dndq] = norm_deriv_dq(xl,yl,L,abar)

  d2wdxdq = abar*(1-cos(2*pi*yl/L))*((2*pi/L)*sin(2*pi*xl/L));
  d2wdydq = abar*(1-cos(2*pi*xl/L))*((2*pi/L)*sin(2*pi*yl/L));

  dndq = [-d2wdxdq; -d2wdydq; 0];

  %dwdx = @(x,y) abar .* ( (cos(2*pi.*y./L)-1).*( (-2*pi./L).*sin(2*pi.*x./L) ) );
  %dwdy = @(x,y) abar .* ( (cos(2*pi.*x./L)-1).*( (-2*pi./L).*sin(2*pi.*y./L) ) );
  %n1_anal = @(x,y) 0.*(x.^0).*(y.^0);
  %n2_anal = @(x,y) 0.*(x.^0).*(y.^0);
  %n3_anal = @(x,y) 1*(x.^0).*(y.^0);
  %dn1dq_anal = @(x,y) -1.*dwdx(x,y);
  %dn2dq_anal = @(x,y) -1.*dwdy(x,y);
  %dn3dq_anal = @(x,y) 0.*(x.^0).*(y.^0);


end
