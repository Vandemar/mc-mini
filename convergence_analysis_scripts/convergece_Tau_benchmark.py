from math import log, sqrt, exp, cos, sin, pi
import h5py
import itertools
import numpy

data_groups = [
 "Divergence",
 "Pressure",
 "Temperature",
 "UVelocity",
 "VVelocity",
 "Viscosity"
 ]

#def v_x_exact(i,j,h):
 #return cos(i*h_x)*sin(j*h_y) 

#def v_y_exact(i,j):
 #return -sin(i*h_x)*cos(j*h_y) 

#def p_exact(i,j,h):
 #return sin(i*h_x)*sin(j*h_y) 

def test():

 ## Path to the HDF5 database in which the computational data is stored.   
 output_data_directory = "../output/tauBenchmark/"

 dataset = []
 dataset_path = []
 data = [] 

 for k in range(0,7):
  dataset_path.append(output_data_directory + "tauBenchmark0" + str(k+1) + "-1.h5")  
  dataset.append(h5py.File(dataset_path[k], "r"))
  data.append(numpy.array(dataset[k]["UVelocity"]))
  
  (L_x, L_y) = (2*pi, pi)
  
  dimensions = (M, N) = data[k].shape
  (h_x, h_y) = (L_x/M, L_y/N)
  
  exact_solution = numpy.empty(dimensions)
  
  for i,j in itertools.product(range(0, M), range(0,N)):  
   exact_solution[i,j] = sin(i*h_x)*sin(j*h_y)

  error.append(data[k] - exact_solution)
  
  if k == 6:
   print("Computed solution = " + str(data[k][0,128]))
   print("Exact solution = " + str(exact_solution[0,0]))
   print(data[k][0,128] - exact_solution[0,0])
   print("Computed solution max = " + str(grid_norm(data[k], numpy.inf)))
  #print(abs(data[k-1]).max()) 
  #print(h_x/h_y) 

 for k in range(0,6):
  error_1 = grid_norm(error[k], 2)
  error_2 = grid_norm(error[k+1], 2)
  #print(error_1)
  #print(error_2)
  #print(error_1/error_2)
  print((log(error_1) - log(error_2))/log(2)) 


def grid_norm(x, norm_type):
 """!
 @brief Grid norm of a function, represented as a vector in @f$ \mathbf{R}^M @f$ (redefined based on <a href="https://github.com/numpy/numpy/blob/v1.9.1/numpy/linalg/linalg.py#L1924">linalg.py</a>)    
 
 @param x Input array.

 @param norm_type Type of the norm, which can be {1, 2, numpy.inf}, i.e. one-, two-, and sup-norm.

 @return Norm of the vector.
 """

# Treat "x" as a one-dimensional numpy array, on which standard accumulation and maximum/minimum search functions can be used.    
 x = numpy.ravel(x)
 h = 1/x.size
    
 if norm_type == numpy.inf:
  return abs(x).max()
 elif norm_type == 1:
  return numpy.sum(abs(x))*h
 elif norm_type == 2:
  return numpy.sqrt(numpy.dot(x, x)*h)