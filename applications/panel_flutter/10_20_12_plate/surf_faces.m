
close all;
clear all;
clc;

% Import elements and nodes
elems = importdata('./elements.txt');
elems = elems.data;
nodes = dlmread('./nodes.txt');
surf_nodes = dlmread('./face_nodes.txt');

surf_node_id = surf_nodes(:,1);
num_elems = max(size(elems));

face_nodes = zeros(num_elems, 11);

parfor (k = 1:num_elems, 12)
  loc_nodes = elems(k,:);
  for j = 1:10
    loc_node = loc_nodes(j);
    onSurf = find(surf_nodes == loc_node);
    if isempty(onSurf)
      loc_nodes(j) = 0;
    end
  end
  if (nnz(loc_nodes) == 6)
    face_nodes(k,:) = [k, loc_nodes];
  end
end

dlmwrite('./surf_faces.txt', face_nodes, 'precision', 8);
