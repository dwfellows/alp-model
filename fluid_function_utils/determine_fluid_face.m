function [fluid_face, struct_bary] = determine_fluid_face(fluid_face_conn, fluid_node_locs, fluid_nodal_faces, structural_node_loc)

  % Find closest fluid node
  vec_diff = structural_node_loc' - fluid_node_locs(:,2:4)';
  vec_diff_mag = vecnorm(vec_diff);

  % Pick off faces connected to this node
  [~,min_node_inds] = sort(vec_diff_mag);

  for j = 1:length(min_node_inds)
    min_node_ind = min_node_inds(j);
    min_node_id = fluid_node_locs(min_node_ind,1);
    min_node_faces = fluid_nodal_faces(min_node_id+1,:);
    min_node_faces = min_node_faces(min_node_faces~=0);
    for i = 1:length(min_node_faces)
      local_fluid_face = min_node_faces(i);
      local_fluid_face_nodes = fluid_face_conn(local_fluid_face,:) + 1;
      local_fluid_face_nodelocs = fluid_node_locs(local_fluid_face_nodes,2:4);

      % Compute projection of point onto plane - for complex geometries
      v1 = local_fluid_face_nodelocs(2,:) - local_fluid_face_nodelocs(1,:);
      v2 = local_fluid_face_nodelocs(3,:) - local_fluid_face_nodelocs(1,:);
      loc_normal = cross(v1,v2);
      structural_node_loc_proj = structural_node_loc - dot(structural_node_loc - local_fluid_face_nodelocs(1,:), loc_normal).*loc_normal;

      %figure; hold on;
      %for k = 1:3
      %  scatter3(local_fluid_face_nodelocs(k,1), local_fluid_face_nodelocs(k,2), local_fluid_face_nodelocs(k,3), 'r');
      %end
      %scatter3(structural_node_loc(1), structural_node_loc(2), structural_node_loc(3), 'b'); pause;

      local_struc_barycentric = fluid_cart2bary(local_fluid_face_nodelocs, structural_node_loc_proj);
      isInside = 1;
      for k = 1:3
        if (local_struc_barycentric(k) < -1e-9)
          isInside = 0;
        end
      end
      if (isInside == 1)
        fluid_face = local_fluid_face;
        struct_bary = local_struc_barycentric;
        return;
      end
    end

  end

end
