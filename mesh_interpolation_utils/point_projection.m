function [point_proj] = point_projection(point_loc, cell_point_locs)
  % Obtain projection of point_loc onto plane defined by cell_point_locs
  n = plane_normal(cell_point_locs);
  p1 = cell_point_locs(1,:);
  point_proj = point_loc - dot(n, point_loc-p1).*n;

end

