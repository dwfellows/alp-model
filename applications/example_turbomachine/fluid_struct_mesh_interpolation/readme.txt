This driver function relies on structural mesh and fluid solution/mesh data obtained from the ANSYS software.

The structural data is derived in general by:
  hp_elems.txt -- export element data from ANSYS meshing software.

  hp_node_elem.txt -- derived from element and node data (obtained from ANSYS meshing software) 
                      by subroutines in the utilities.

  hp_nodes.txt -- export node data from ANSYS meshing software.

  hp_face_nodes.txt -- export node data on compliant faces only from ANSYS meshing software.

  hp_faces.txt -- derived from node data, element data, and face_node data obtained from ANSYS
                    meshing software utilizing subroutines in the utilities.

  mesh_ref_pts_hpt.csv -- reference points of structural mesh used to align structural and fluid meshes.
                            Selected by hand from ANSYS meshing software.

The fluid data is derived in general by:
  hp_face_data.csv -- Export flow solution quantities from CFD-Post, and delete point solution data 
                       (keep only face/point connectivity).

  hp_pt_data.csv -- Export flow solution quantities from CFD-Post, and delete point/face data
                      (keep only point solution data).
