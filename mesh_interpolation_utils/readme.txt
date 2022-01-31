Mesh-to-mesh interpolation utilities -- these were originally developed separately from the other utilities making up the ALP model.
These are written for utilization with HP-turbine data (as developed on), may need to be modified for LP-turbine data.

As a consequence, there may be some discrepancies/bugs that have arisen when modifying the source code to "fit in" to a utilities folder.

The main routines are:
  face_to_pt.m -- Take in ANSYS face data and construct face/pt connectivity needed for mesh-to-mesh interpolation routines.

  construct_integration_locs.m -- Take in ANSYS structural mesh data and derive locations of integraiton points for spatial integration
                                   on structural mesh.

  fluid_mesh_interpolation.m -- Take in strucctural mesh data and fluid solution data and interpolation data from fluid mesh
                                 onto structural integration points. Writes file of fluid data and integration point # for use
                                 in compliant surface integration routines.

All other routines are helper functions for the above main routines.
