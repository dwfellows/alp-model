function [bool_contain] = triangle_test(min_node_loc, neighbor1_node_loc, neighbor2_node_loc, fluid_node_loc)

  full_area = triangle_area(min_node_loc, neighbor1_node_loc, neighbor2_node_loc);

  area1 = triangle_area(fluid_node_loc, min_node_loc, neighbor1_node_loc);
  area2 = triangle_area(fluid_node_loc, min_node_loc, neighbor2_node_loc);
  area3 = triangle_area(fluid_node_loc, neighbor1_node_loc, neighbor2_node_loc);

  if ((area1 + area2 + area3) == full_area)
    bool_contain = 1;
  else
    bool_contain = 0;
  end

%  % Construct vectors
%  min_node_vec = min_node_loc - fluid_node_loc;
%  neighbor1_node_vec = neighbor1_node_loc - fluid_node_loc;
%  neighbor2_node_vec = neighbor2_node_loc - fluid_node_loc;
%
%  % Find angles
%  min_node_neighbor1_cosang = max(min(dot(min_node_vec,neighbor1_node_vec)/(norm(min_node_vec)*norm(neighbor1_node_vec)),1),-1);
%  min_node_neighbor1_angd = real(acosd(min_node_neighbor1_cosang));
%
%  min_node_neighbor2_cosang = max(min(dot(min_node_vec,neighbor2_node_vec)/(norm(min_node_vec)*norm(neighbor2_node_vec)),1),-1);
%  min_node_neighbor2_angd = real(acosd(min_node_neighbor2_cosang));
%
%  n1_n2_cosang = max(min(dot(neighbor1_node_vec,neighbor2_node_vec)/(norm(neighbor1_node_vec)*norm(neighbor2_node_vec)),1),-1);
%  n1_n2_angd = real(acosd(n1_n2_cosang));
%
%  combined_angle = min_node_neighbor1_angd + min_node_neighbor2_angd + n1_n2_angd;
%  triple_prod = dot(min_node_vec, cross(neighbor1_node_vec,neighbor2_node_vec));
%
%
%  if ((combined_angle > 359.9) & (triple_prod == 0))
%    bool_contain = 1;
%  else
%    bool_contain = 0;
%  end

end
