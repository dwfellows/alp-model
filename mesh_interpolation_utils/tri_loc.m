% Triangular interpolation function for location
function [loc] = tri_loc(face_node_locs, interp_loc)
  % Linear triangle location from area coordinates
  % face_node_locs: Cartesian locations of each face node
  % interp_loc: location of point in triangle, area coordinates
  % loc: location of point in triangle, Cartesian coordinates

  loc = face_node_locs(1,:).*interp_loc(1) + face_node_locs(2,:).*interp_loc(2) + face_node_locs(3,:).*interp_loc(3);

end

