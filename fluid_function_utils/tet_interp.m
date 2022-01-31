function [interp_val] = tet_interp(node_vals, bary_coords)

  L1 = bary_coords(1);
  L2 = bary_coords(2);
  L3 = bary_coords(3);
  L4 = bary_coords(4);

  interp_val = node_vals(1,:).*((2*L1-1)*L1) + node_vals(2,:).*((2*L2-1)*L2) ...
             + node_vals(3,:).*((2*L3-1)*L3) + node_vals(4,:).*((2*L4-1)*L4) ...
             + node_vals(5,:).*(4*L1*L2) + node_vals(6,:).*(4*L2*L3) ...
             + node_vals(7,:).*(4*L1*L3) + node_vals(8,:).*(4*L1*L4) ...
             + node_vals(9,:).*(4*L2*L4) + node_vals(10,:).*(4*L3*L4);

end
