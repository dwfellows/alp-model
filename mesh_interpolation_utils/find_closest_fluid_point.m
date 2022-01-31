function [closest_fluid_point] = find_closest_fluid_point(struct_point, surf_points)
  % Returns fluid point in surf_poitns that is closest to struct_point

  vect_dist = surf_points - struct_point;
  mag_dist = vecnorm(vect_dist, 2, 2);
  [min_val, min_ind] = min(mag_dist);
  closest_fluid_point = min_ind;

end

