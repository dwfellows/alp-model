function [R] = rotation_matrix_align(a,b)
  % Rotation matrix which aligns a onto b

  v = cross(a,b); c = dot(a,b);
  if (norm(v) == 0)
    R = eye(3);
  else
    vx = [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];
    R = eye(3) + vx + vx*vx.*(1/(1+c));
  end

end

