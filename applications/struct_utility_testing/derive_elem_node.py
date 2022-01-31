
import numpy as np
import scipy
import pyansys

def construct_elem_node_dict(rst_filename):

  # Read binary file
  binary_file = pyansys.read_binary(rst_filename)

  # Obtain element file
  elem_nodes = binary_file.geometry.elem
  elem_nodes = np.array(elem_nodes)

  # Strip unwanted data -- this part only keeps 10 nodes, quadratic tet
  elem_nodes = elem_nodes[:,10:20]

  # Write to csv file
  csv_name = rst_filename[0:-4] + '.csv'
  np.savetxt(csv_name, elem_nodes, delimiter=',')
