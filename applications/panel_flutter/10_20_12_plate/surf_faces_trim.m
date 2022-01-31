
close all;
clear all;
clc;

% Import elements and nodes
faces = dlmread('./surf_faces.txt');

num_elems = max(size(faces));
face_nodes = [];
for k = 1:num_elems
  disp(k);
  if (nnz(faces(k,:)) > 0)
    face_nodes = [face_nodes; faces(k,:)];
  end
end

dlmwrite('./surf_faces_trim.txt', face_nodes, 'precision', 8);
