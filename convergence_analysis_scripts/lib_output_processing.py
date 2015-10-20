"""!
  @namespace lib_output_processing
  
  Python scripts devoted to estimating convergence rates of the numerical solutions.

  Convergence is estimated through two methods:
  1. Direct comparison against the known exact solution for the Tau benchmark problem (cf. Tau, Numerical Solution of the Steady Stokes Equations).
  TODO: 2. Richardson's extrapolation (NOTE: What is the proper name for this algorithm that would be recognized by most reserachers in the field?).

  TODO: Redesign this by using classes since the data and functions required for the two methods largely overlap.
  
"""

from math import log, sqrt, cos, sin, pi
import h5py
import itertools
import numpy

# All possible datasets provided by mc-mini. 
datasets = [
 "Divergence",
 "Pressure",
 "Temperature",
 "UVelocity",
 "VVelocity",
 "Viscosity"
 ]

# Datasets relevant to Tau benchmark computation.
tau_benchmark_datasets = [
 "Pressure",
 "UVelocity",
 "VVelocity"
 ]

# Header strings for the convergence table.
convergence_table_header_data = { 
 "k" : "k",
 "M" : "M = 2*2^k",
 "N" : "N = 2^k", 
 "h" : "h = 2*pi/M", 
 "sup_norm" : "Sup-Norm", 
 "sup_norm_rate" : "Rate", 
 "one_norm" : "One-Norm", 
 "one_norm_rate" : "Rate", 
 "two_norm" : "Two-Norm", 
 "two_norm_rate" : "Rate"
 }

# Location of the output hdf5 files.
output_data_directory = "../output"

def tau_benchmark_convergence():
 """!
 @brief Compute convergence rates for all Tau benchmark datasets in various norms and output them into a table.
 
 @return Outputs a table containing convergence rates' estimates.
 """
 
 # Path to the text file where the table is to be stored.
 table_file_path = "." + "/convergence_tables/test.txt"
 
 # Write the table to a text file.  
 output_text_table_file = open(table_file_path, "w") 
 
 for dataset in tau_benchmark_datasets:
  
  # Compute norms of the errors.
  error_norms = error_norms_data(dataset)
 
  # Formatted string containing the attributes of the computation.
  opening = "Dataset: {dataset_name}".format(
   dataset_name = dataset
   )
  
  #Formatted string containing the text of the header of the table.
  header = "{k:^3} {M:^10} {N:^10} {h:^15} {sup_norm:^15} {sup_norm_rate:^5} {one_norm:^15} {one_norm_rate:^5} {two_norm:^15} {two_norm_rate:^5}".format(**convergence_table_header_data)
  
  # Return a formatted string containing convergence information for the given grid resolution. 
  def convergence_row(k): 
   return "{k:^3} {M:^10} {N:^10} {h:^15.10f} {sup_norm:^15.10f} {sup_norm_rate:^5.2f} {one_norm:^15.10f} {one_norm_rate:^5.2f} {two_norm:^15.10f} {two_norm_rate:^5.2f}".format(**convergence_table_row_data(k, error_norms))

  # Write the table to the text file.  
  output_text_table_file.write(opening + "\n\n")
  output_text_table_file.write(header + "\n\n")
  for k in range(1, 3):
   output_text_table_file.write(convergence_row(k) + "\n")
  output_text_table_file.write("\n")
 
 return None

def grid_norm(x, norm_type):
 """!
 @brief Grid norm of a function on a discrete grid, represented as a vector in @f$ \mathbf{R}^M @f$ (redefined based on <a href="https://github.com/numpy/numpy/blob/v1.9.1/numpy/linalg/linalg.py#L1924">linalg.py</a>)    
 
 @param x Input array.

 @param norm_type Type of the norm, which can be {1, 2, numpy.inf}, i.e. one-, two-, and sup-norm.

 @return Norm of the vector.
 """

# Treat "x" as a one-dimensional numpy array, on which standard accumulation and maximum/minimum search functions can be used.    
 x = numpy.ravel(x)
 h = 2*pi*pi/x.size
    
 if norm_type == numpy.inf:
  return abs(x).max()
 elif norm_type == 1:
  return numpy.sum(abs(x))*h
 elif norm_type == 2:
  return numpy.sqrt(numpy.dot(x, x)*h)

def error_norms_data(dataset): 
 """!
 @brief Compute and store various norms of the error vector for all grid resolutions for a given dataset.
 
 @param dataset One of the valid datasets. @see datasets
 
 @return Array of dictionaries indexed by grid refinement exponent containing "norm : error" key-value pairs
 """
 
 # Domain dimensions.
 (L_x, L_y) = (2*pi, pi)

# Exact solution for pressure.  
 def p_exact(i,j):
  return sin( (i+0.5)*h_x )*sin( (j+0.5)*h_y ) 

# Exact solution for uVelocity restricted to the cell-centers of the grid by interpolation (averaging).
 def v_x_exact_interpolated(i,j):
   # NOTE: This is the correct index ordering, i.e. "ij" <=> "xy"
  return 0.5 * ( cos(i*h_x)*sin( (j+0.5)*h_y) + cos( (i+1)*h_x )*sin( (j+0.5)*h_y))
  # NOTE: This is the incorrect index ordering, i.e. "ji" <=> "xy"
  # Convergence takes place for this (incorrect) index ordering.
  #return 0.5 * ( cos((j+0.5)*h_y)*sin( i*h_x) + cos( (j+0.5)*h_y )*sin((i+1)*h_x ))

# Exact solution for vVelocity restricted to the cell-centers of the grid by interpolation (averaging).
 def v_y_exact_interpolated(i,j):
  # NOTE: This is the correct index ordering, i.e. "ij" <=> "xy"
  return 0.5 * (-sin( (i+0.5)*h_x )*cos(j*h_y) - sin( (i+0.5)*h_x )*cos( (j+1)*h_y))
  # NOTE: This is the incorrect index ordering, i.e. "ji" <=> "xy"
  # Convergence takes place for this (incorrect) index ordering.
  #return 0.5 * (-sin( j*h_y )*cos((i+0.5)*h_x) - sin( (j+1)*h_y )*cos( (i+0.5)*h_x)) 

# Dictionary containing "dataset_name : corresponding_exact_solution" key-value pairs.
 exact_solutions = {
  "Pressure"  : p_exact,
  "UVelocity" : v_x_exact_interpolated,
  "VVelocity" : v_y_exact_interpolated
 }
 
 # Path to the directory where the computational data is stored.   
 output_directory = output_data_directory + "/tauBenchmark/"

 files = []
 file_path = []
 data = [] 
 errors = []

# Append files to an array
 for k in range(0,3):  
  file_path.append(output_directory + "tauBenchmark0" + str(k+1) + "-1.h5")  
  files.append(h5py.File(file_path[k], "r"))
 
 # Array of dictionaries containing errors in various norms
 _error_norms = {}

 for k in range(0,3):  
  data.append(numpy.array(files[k][dataset]))

  dimensions = (M, N) = data[k].shape
  (h_x, h_y) = (L_x/N, L_y/M)
  
  exact_solution = numpy.empty(dimensions)
  for i,j in itertools.product(range(0, M), range(0, N)):  
   exact_solution[i,j] = exact_solutions[dataset](i,j)
  
  # Error vector
  error = (data[k] - exact_solution)
  # Store various norms of the error vector  
  _error_norms[k] = {
  "sup_norm" : grid_norm(error, numpy.inf), 
  "one_norm" : grid_norm(error, 1), 
  "two_norm" : grid_norm(error, 2)
  }
 
 # Close the files 
 for k in range(0,3):    
  files[k].close() 

 return _error_norms

def convergence_table_row_data(k, _error_norms):
 """!
 @brief Compute and store various norms of the error vector for a given grid resolution.

 @param k Refinement exponent in the expression for grid stepsize: h = 2*pi/2*2**(k+2)
 
 @param _error_norms Array of dictionaries indexed by grid refinement exponent containing "norm : error" key-value pairs. 
 
 @return Convergence data in the form of a dictionary containing key-value pairs necessary to form one row in the convergence table 
 for the given refinement exponent.
 """ 
 
# Convergence rate for the computation is an estimate of the order of convergence of the numerical method, 
# and it is computed as log_2(error[twice smaller grid resolution]/error[current grid resolution])
 def convergence_rate (k, norm): 
  return abs(log(_error_norms[k-1][norm]) - log(_error_norms[k][norm]))/log(2)

# Store the relevant data in the dictionary 
 data = {}
 data["k"] = k+2
 data["M"] = 2*2**(k+2)
 data["N"] = 2**(k+2) 
 data["h"] = 2*pi/data["M"] 
 data["sup_norm"] = _error_norms[k]["sup_norm"] 
 data["sup_norm_rate"] = convergence_rate(k, "sup_norm") 
 data["one_norm"] = _error_norms[k]["one_norm"]
 data["one_norm_rate"] = convergence_rate(k, "one_norm")
 data["two_norm"] = _error_norms[k]["two_norm"]
 data["two_norm_rate"] = convergence_rate(k, "two_norm")
 return data