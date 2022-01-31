function [L1, L2, L3] = barycentric(point_loc, cell_point_locs)

  n = plane_normal(cell_point_locs);
  e01 = cell_point_locs(2,:) - cell_point_locs(1,:);
  e12 = cell_point_locs(3,:) - cell_point_locs(2,:);
  e20 = cell_point_locs(1,:) - cell_point_locs(3,:);
  A = 0.5*dot(cross(e01,e12), n);
  Au = 0.5*dot(cross(e12, point_loc - cell_point_locs(2,:)), n);
  Av = 0.5*dot(cross(e20, point_loc - cell_point_locs(3,:)), n);

  L1 = Au / A;
  L2 = Av / A;
  L3 = 1 - L1 - L2;

end

