
close all;
clear all;
clc;

elem_node = importdata('./elements.txt');
elem_node = elem_node.data;

num_elems = max(size(elem_node));
num_nodes = max(max(elem_node));
node_elem = zeros(num_nodes, 60);
for k = 1:num_elems
  for j = 1:10
    loc_node = elem_node(k,j);
    i = 0;
    while (i>-1)
      i = i+1;
      test = node_elem(loc_node, i);
      if (test==0)
        node_elem(loc_node, i) = k;
        i = -2;
      end
    end
  end
end

dlmwrite('./node_elem.txt', node_elem, 'precision', 8);
