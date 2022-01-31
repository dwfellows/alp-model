function [] = plot_tet(loc_nodes)

  plot3(loc_nodes(1,1), loc_nodes(1,2), loc_nodes(1,3), loc_nodes(5,1), loc_nodes(5,2), loc_nodes(5,3));
  plot3(loc_nodes(1,1), loc_nodes(1,2), loc_nodes(1,3), loc_nodes(7,1), loc_nodes(7,2), loc_nodes(7,3));
  plot3(loc_nodes(1,1), loc_nodes(1,2), loc_nodes(1,3), loc_nodes(8,1), loc_nodes(8,2), loc_nodes(8,3));
  plot3(loc_nodes(2,1), loc_nodes(2,2), loc_nodes(2,3), loc_nodes(5,1), loc_nodes(5,2), loc_nodes(5,3));
  plot3(loc_nodes(2,1), loc_nodes(2,2), loc_nodes(2,3), loc_nodes(6,1), loc_nodes(6,2), loc_nodes(6,3));
  plot3(loc_nodes(2,1), loc_nodes(2,2), loc_nodes(2,3), loc_nodes(9,1), loc_nodes(9,2), loc_nodes(9,3));
  plot3(loc_nodes(3,1), loc_nodes(3,2), loc_nodes(3,3), loc_nodes(6,1), loc_nodes(6,2), loc_nodes(6,3));
  plot3(loc_nodes(3,1), loc_nodes(3,2), loc_nodes(3,3), loc_nodes(7,1), loc_nodes(7,2), loc_nodes(7,3));
  plot3(loc_nodes(3,1), loc_nodes(3,2), loc_nodes(3,3), loc_nodes(10,1), loc_nodes(10,2), loc_nodes(10,3));
  plot3(loc_nodes(4,1), loc_nodes(4,2), loc_nodes(4,3), loc_nodes(8,1), loc_nodes(8,2), loc_nodes(8,3));
  plot3(loc_nodes(4,1), loc_nodes(4,2), loc_nodes(4,3), loc_nodes(9,1), loc_nodes(9,2), loc_nodes(9,3));
  plot3(loc_nodes(4,1), loc_nodes(4,2), loc_nodes(4,3), loc_nodes(10,1), loc_nodes(10,2), loc_nodes(10,3));

end

