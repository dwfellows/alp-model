function [n] = plane_normal(cell_point_locs)
  % Define the normal to the plane defined by triangle vertices
  p1 = cell_point_locs(1,:);
  p2 = cell_point_locs(2,:);
  p3 = cell_point_locs(3,:);

  vec1 = p2 - p1;
  vec2 = p3 - p1;
  n = cross(vec1, vec2);
  n = n ./ norm(n);

end

