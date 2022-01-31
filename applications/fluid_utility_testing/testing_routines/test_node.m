
close all;
clear all;
clc;

fs = 18;

% Add tools path
addpath('../../../fluid_function_utils');

% Read in test node
tnodes = dlmread('./node_locs.csv');

% Read in circular data
coarse_elem_node = dlmread('../fem_data/circular_coarse/coarse.csv');
coarse_node_elem = dlmread('../fem_data/circular_coarse/coarse_node_elem.csv');
coarse_xdir_deform = dlmread('../fem_data/circular_coarse/coarse_xdir_deformation.txt');
coarse_ydir_deform = dlmread('../fem_data/circular_coarse/coarse_ydir_deformation.txt');
coarse_zdir_deform = dlmread('../fem_data/circular_coarse/coarse_zdir_deformation.txt');
coarse_node_locs = coarse_xdir_deform(:,2:4);
coarse_surface_node_data = dlmread('../fem_data/circular_coarse/coarse_surface.txt');

% Check random nodes and ensure (by eye) the returned element is correct
rand_num = 10;
for l = 1:rand_num
  r = rand(1) * 0.1;
  t = rand(1) * 2*pi;
  x = r*cos(t);
  y = r*sin(t);

  tnode = [x, y, 0.0001];
  [struc_elem] = determine_elem(coarse_elem_node, coarse_node_elem, coarse_surface_node_data, coarse_node_locs, tnode);
  struc_elem_nodes = coarse_elem_node(struc_elem,:);
  struc_elem_locs = coarse_node_locs(struc_elem_nodes,:);
  [fluid_node_loc_barycentric] = cart2bary(struc_elem_locs, tnode);
  disp(vpa(fluid_node_loc_barycentric));

  figure;
  for k = 1:10
    node = struc_elem_nodes(k);
    node_loc = coarse_node_locs(node,:);
    scatter3(node_loc(1), node_loc(2), node_loc(3), 'b'); hold on;
    text(node_loc(1), node_loc(2), node_loc(3), num2str(k)); hold on;
  end
  scatter3(tnode(1), tnode(2), tnode(3), 'r');

  [normal] = determine_normal(struc_elem_locs, tnode, fluid_node_loc_barycentric);
  vpa(normal)
  normal = normal ./ 10000;
  quiver3(tnode(1), tnode(2), tnode(3), normal(1), normal(2), normal(3));
  

  pause;
end

