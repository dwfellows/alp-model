
% Takes ANSYS face connectivity data and inverses the matrix
% to construct an array denoting which faces contain each point

function [message] = face_to_pt(fluid_info)

  message = 0;

  % Read in face data
  face_data_fname = fluid_info(1);
  % face_data_fname = 'hp_face_data.csv';
  face_data = dlmread(face_data_fname);
  
  % Find number of points
  n_pts = max(max(face_data)) + 1;
  n_faces = max(size(face_data));
  
  % Allocate array
  point_faces = zeros(n_pts,10);
  
  % Derive point-face data
  for i = 1:n_faces
    disp(i);
    loc_face = face_data(i,:);
    for j = 1:3
      loc_point = loc_face(j) + 1;
      ent = find(point_faces(loc_point,:)==0, 1, 'first');
      point_faces(loc_point, ent) = i;
    end
  end
  
  dlmwrite('../computation_data/fluid_solution_data/hp_point_face_data.csv', point_faces);
  message = 1;

end
