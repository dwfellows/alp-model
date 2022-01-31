function [barycentric] = fluid_cart2bary(fluid_nodelocs, struct_cart);

  %barycentric = (fluid_nodelocs')\(struct_cart');

  v1 = fluid_nodelocs(2,:) - fluid_nodelocs(1,:);
  v2 = fluid_nodelocs(3,:) - fluid_nodelocs(1,:);
  n = cross(v1,v2);
  a = n ./ norm(n); a = a';
  b = [0;0;1];
  v = cross(a,b);
  c = dot(a,b);
  vx = [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];
  if (c==-1)
    R = eye(3);
  else
    R = eye(3) + vx + vx.*vx.*(1/(1+c));
  end

  fluid_locs = R*(fluid_nodelocs'); fluid_locs = fluid_locs(1:2,:);
  struc_loc = R*struct_cart'; struc_loc = struc_loc(1:2);

  struc_loc = struc_loc - fluid_locs(:,3);
  fluid_locs = [fluid_locs(:,1) - fluid_locs(:,3), fluid_locs(:,2) - fluid_locs(:,3)];

  barycentric = fluid_locs \ struc_loc;
  barycentric = [barycentric; 1 - barycentric(1) - barycentric(2)];

end
