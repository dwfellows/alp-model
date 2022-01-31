function [area] = triangle_area(p1, p2, p3)

  vec1 = p2 - p1;
  vec2 = p3 - p1;

  vec3 = cross(vec1, vec2);
  area = norm(vec3) / 2;

end
