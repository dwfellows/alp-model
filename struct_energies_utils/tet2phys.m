function [physical_loc] = tet2phys(tet_loc, x1, x2, x3, x4)

  L1 = tet_loc(1);
  L2 = tet_loc(2);
  L3 = tet_loc(3);
  L4 = tet_loc(4);

  x = L1*x1(1) + L2*x2(1) + L3*x3(1) + L4*x4(1);
  y = L1*x1(2) + L2*x2(2) + L3*x3(2) + L4*x4(2);
  z = L1*x1(3) + L2*x2(3) + L3*x3(3) + L4*x4(3);

  physical_loc = [x,y,z];

end

