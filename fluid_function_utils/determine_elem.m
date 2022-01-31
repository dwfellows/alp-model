function [elem] = determine_elem(elem_node, node_elem, surface_node_data, structural_node_locs, fluid_node_loc)
  % Assumes that all node locations (between structure and fluid meshes) are with respect to the same coordinate axis

  % Surface nodes: 1 NodeID | 2 X-loc | 3 Y-loc | 4 Z-loc
  surface_node_locs = surface_node_data(:,2:4);
  surface_node_ids = surface_node_data(:,1);

  % Find closest surface nodes
  vec_diff = fluid_node_loc' - surface_node_locs';
  vec_diff_mag = vecnorm(vec_diff);

  % Pick off elements connected to this node
  [~,min_node_inds] = sort(vec_diff_mag);
  for j = 1:length(min_node_inds)
    min_node_ind = min_node_inds(j);
    min_node_id = surface_node_data(min_node_ind,1);
    min_node_loc = structural_node_locs(min_node_ind,:);

    min_node_contained_elems = node_elem(min_node_id,:);
    min_node_contained_elems = min_node_contained_elems(min_node_contained_elems~=0);
  
    % Check for containment on surface
    for i = 1:length(min_node_contained_elems)
      % Obtain all nodes in particular element
      loc_elem = min_node_contained_elems(i);
      loc_elem_nodes = elem_node(loc_elem,:);

      % Get nodal locations
      loc_elem_nodes_locs = structural_node_locs(loc_elem_nodes, :);
      
      % Acquire barycentric coordinates of fluid node location
      [fluid_node_loc_barycentric] = cart2bary(loc_elem_nodes_locs, fluid_node_loc);
%      disp(fluid_node_loc_barycentric);
%      figure(1); clf;
%      scatter3(loc_elem_nodes_locs(:,1), loc_elem_nodes_locs(:,2), loc_elem_nodes_locs(:,3), 'b'); hold on;
%      scatter3(fluid_node_loc(1), fluid_node_loc(2), fluid_node_loc(3)); 
%      pause;

      % Check barycentric coordinates and determine element containment
      isInside = 1;
      for k = 1:4
        if (fluid_node_loc_barycentric(k) < -1e-9)
          isInside = 0;
        end
      end

      if (isInside == 1)
        elem = loc_elem;
        return;
      end

    end
  end

end

