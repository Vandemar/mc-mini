"""!
  @namespace lib_output_processing
  
  Python scripts devoted to processing the computational results, such as:
  <ul>
    <li>
      creating the computational database and preparing the appropriate computational context
    </li>
    <li>
      calculating the convergence rates
    </li>
    <li>
      plotting the results
    </li> 
    <li>
      incorporating the processed data into HTML documents for greater readability
    </li>
  </ul>
"""

import os.path

from math import log, sqrt, exp, cos, sin, pi

import h5py

import itertools

import matplotlib
import matplotlib.pyplot as pyplot

import numpy
from numpy.core import asarray

from lxml import etree

#class convergence(object):
 #h =   

data_groups = [
 "Divergence",
 "Pressure",
 "Temperature",
 "UVelocity",
 "VVelocity",
 "Viscosity"
 ]

def v_x_exact(i,j,h):
 return cos(i*h)*sin(j*h) 

def v_y_exact(i,j,h):
 return -sin(i*h)*cos(j*h) 

def p_exact(i,j,h):
 return sin(i*h)*sin(j*h) 

def test():

 ## Path to the HDF5 database in which the computational data is stored.   
 output_data_directory = "../output/tauBenchmark/"

 files = []
 dataset = []
 dataset_path = []
 data = [] 

 error = []
  
 h_x = []
 h_y = []
 
 for k in range(0,7):
  dataset_path.append(output_data_directory + "tauBenchmark0" + str(k+1) + "-1.h5")  
  dataset.append(h5py.File(dataset_path[k], "r"))
  data.append(numpy.array(dataset[k]["UVelocity"]))
  
  (L_x, L_y) = (2*pi, pi)
  
  ## NOTE: M and N are flipped in the code.
  #(_N, _M) = data[k].shape
  dimensions = (M, N) = data[k].shape
  (h_x, h_y) = (L_x/N, L_y/M)
  
  exact_solution = numpy.empty(dimensions)
  
  for i,j in itertools.product(range(0, M), range(0, N)):  
   exact_solution[i,j] = cos((i+0.5))*sin((j+0.5)*h_y)
  
  #print(h_x)
  #print(h_y)
  #print(M)
  #print(N)
  
  error.append(data[k] - exact_solution)
  
  #if k == 6:
   #print("Computed solution = " + str(data[k][0,128]))
   #print("Exact solution = " + str(exact_solution[0,0]))
   #print(data[k][0,128] - exact_solution[0,0])
   #print("Computed solution max = " + str(grid_norm(data[k], numpy.inf)))
  ##print(abs(data[k-1]).max()) 
  ##print(h_x/h_y) 

 for k in range(0,6):
  error_1 = grid_norm(error[k], 2)
  error_2 = grid_norm(error[k+1], 2)
  #print(error_1)
  #print(error_2)
  #print(error_1/error_2)
  print((log(error_1) - log(error_2))/log(2)) 

#h_x = 2*pi/M 
#h_y = pi/N



## Path to the website directory where the generated HTML documents will be written
html_path = "website/"


## Valid initial condition input strings 
initial_conditions = [
 "Square_Wave", 
 "Semicircle", 
 "Gaussian_Pulse"
  ]   
    
## Valid fluxes input strings    
fluxes = [
 "Upwind", 
 "Lax_Friedrichs", 
 "Lax_Wendroff", 
 "Fromm", 
 "Fromm_CFL_half", 
 "Fromm_van_Leer", 
 "Fromm_van_Leer_CFL_half",
 "Flux_Corrected_Transport", 
 "Lax_Wendroff_Fourth_Order"    
  ]

## Header strings for the convergence table
convergence_table_header_data = { 
 "k" : "k",
 "M" : "M = 2^k", 
 "h" : "h = M^(-1)", 
 "sup_norm" : "Sup-Norm", 
 "sup_norm_rate" : "Rate", 
 "one_norm" : "One-Norm", 
 "one_norm_rate" : "Rate", 
 "two_norm" : "Two-Norm", 
 "two_norm_rate" : "Rate"
 }

def initial_data(initial_condition, k):
 """!
 @brief Generate the initial data based on the desired initial condition.
 
 @param initial_condition One of the valid initial conditions. @see initial_conditions

 @param k Refinement exponent in the expression for grid stepsize: h = 2^(-k)

 @return Array containing requested initial data.
 """

 data = []
 
 M = 2**k
 h = 1/M
 
 def square_wave(x):
  if (abs(x - 0.5) <= 0.25):
   return 1
  else:
   return 0
 
 def semicircle(x):
  return sqrt(0.25-pow(x - 0.5, 2))
 
 def gaussian_pulse(x):
  return exp(-256*pow(x - 0.5, 2)) 
 
 if initial_condition == "Square_Wave":
  for i in range(0, M):
   data.append(square_wave(h*i))
 elif initial_condition == "Semicircle":
  for i in range(0, M):
   data.append(semicircle(h*i))
 else:
  for i in range(0, M):
   data.append(gaussian_pulse(h*i))
   
 return data

def create_output_database():
 """! 

 @brief Create an hdf5 file in which the results of the computation will be stored.
    
 Creates an hdf5 file with the necessary group hierarchy to store the results of 
 computations in an organized manner. 
    
 This function should only be used once: once the hdf5 file is created with the 
 appropriate group structure, the main program shall not alter that structure.

  """

# Do not alter the hdf5 file if it already exists
 if os.path.exists(database_path):
  print("DATABASE STATUS:")
  print("\t" + database_path + " already exists and is ready to store the results of computations")
  return None
# Create hdf5 file. The flag "-w" means "create file, fail if exists" 
 else:
  computations_database = h5py.File(database_path, "w-")

# Create initial data datasets and write initial data into them 
  for initial_condition in initial_conditions:
   for k in range (6,17):
    dataset_initial_path = initial_condition + "/k = " + str(k) + " initial_data"
    computations_database[dataset_initial_path] = initial_data(initial_condition, k)
# Create data groups for storing the results of computations 
   for flux in fluxes:     
    group_path = initial_condition + "/" + flux
    computations_database.create_group(group_path)

# Write the appropriate attributes that are needed for particular computations, 
# i.e. create the appropriate environment for each computational method     
    computations_database[group_path].attrs["a"] = 3.0
    computations_database[group_path].attrs["T"] = 9.0
    if flux == "Lax_Wendroff_Fourth_Order": 
     computations_database[group_path].attrs["CFL"] = 0.2
    elif flux in ["Fromm_CFL_0.5", "Fromm_van_Leer_CFL_0.5"]:
     computations_database[group_path].attrs["CFL"] = 0.5
    else:
     computations_database[group_path].attrs["CFL"] = 0.9
     
  computations_database.close()   
  print("DATABASE STATUS:")
  print("\t" + database_path + " has been created and is ready to store the results of computations")
   
def create_plots(initial_condition, flux): 
 """! 

 @brief Plot initial data and four selected computations to illustrate behavior of the numerical method as the grid resolution
 increases.

 @param initial_condition One of the valid initial conditions. @see initial_conditions

 @param flux One of the valid fluxes. @see fluxes

 """
 
# Path to the hdf5 file containing results of computations
 computations_database = h5py.File(database_path, "r")
 
# Path to the data group containing the computational data
 group_path = initial_condition + "/" + flux
 
# Path to the highest grid resolution data (k = 16)
 dataset_path = group_path + "/k = " + str(16)
# Adjust the plot window to the maximum and minimum values the computed solution attains 
 ymax = asarray(computations_database[dataset_path]).max()
 ymin = asarray(computations_database[dataset_path]).min()
 pyplot.ylim(ymax = ymax*1.2)
 pyplot.ylim(ymin = ymin*1.2)
 
# Path to the initial condition dataset 
 dataset_initial_path = initial_condition + "/k = " + str(16) + " initial_data"
 
# The title of the plot, containing relevant information about the computation 
 title = "{group_path}: CFL = {CFL}, a = {a}".format(
 group_path = group_path.replace("/", ", "), 
 CFL = computations_database[group_path].attrs["CFL"], 
 a = computations_database[group_path].attrs["a"]
 )
 
# Define the title and axes' labels; enable grid in the plot's background
 pyplot.title(title)
 pyplot.xlabel("x")
 pyplot.ylabel("u(x, t = 9)")
 pyplot.grid(True)

# Define the grid against which data should be ploted, i.e. uniform grid with stepsize 2^(-16)
 x = numpy.linspace(0, 1, 2**16)
# Plot the initial data 
 pyplot.plot(x, computations_database[dataset_initial_path])

# Include the initial data legend in the array of all legends which will be displayed in the plot
 legend_labels = ["Initial Data"]
# Plot the legend for the initial data 
 pyplot.legend(legend_labels)

# Plot solutions corresponding to grid stepsizes 2^(-10), 2^(-12), 2^(-14), and 2^(-16),
 for k in range(10, 17, 2):
  x = numpy.linspace(0, 1, 2**k)   
  
  dataset_path = group_path + "/k = " + str(k)
 
  pyplot.plot(x, computations_database[dataset_path])
 
# Append the legend text corresponding to the current computation to the array of legends
  legend_labels.append("k = " + str(k))
# Plot the newly appended legend; ncol - number of columns, fancybox - rounded corners of the box enclosing the legend text
  pyplot.legend(legend_labels, ncol = 3, fancybox = True, fontsize = "small")

# Output plots in two formats to be used in the website 
 pyplot.savefig(html_path + "images/" + group_path.replace("/","_") + ".png")
 pyplot.savefig(html_path + "pdf/plots/" + group_path.replace("/","_") + ".pdf")
 pyplot.close()
 
 computations_database.close()
 
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

def check_conservation(initial_condition, flux):
 """!
 @brief Check whether the Riemman intergral of the solution is equal to the integral of the initial data, i.e. check if the numerical
 method preserved the conservation property of the linear advection, up to the floating point error.
 
 @param initial_condition One of the valid initial conditions. @see initial_conditions

 @param flux One of the valid fluxes. @see fluxes
 """
 
# Path to the hdf5 file containing results of computations
 computations_database = h5py.File(database_path, "r")
 
# Path to the data group containing the computational data
 group_path = initial_condition + "/" + flux
 
 for k in range (6, 17):
  dataset_initial_path = initial_condition + "/k = " + str(k) + " initial_data"
  dataset_path = group_path + "/k = " + str(k)
  
  h = 2**(-k)
  initial_data = asarray(computations_database[dataset_initial_path])
  computed_solution = asarray(computations_database[dataset_path])
# Conservation error is the absolute value of the difference between the integrals of the initial data and the computed solution
  conservation_error = abs(numpy.sum(initial_data) - numpy.sum(computed_solution))*h

# Criterion of preserving conservation must take into account the floating point error accumulation during the computation
# Experience shows that 0.1 is the upper bound for the floating point error for a conservative numerical method.
  if (conservation_error > 1e-1):
   print("ERROR: \n\t" + dataset_path + ": conservation is not preserved")
   print("\tConservation floating-point error: {0:.1e}".format(conservation_error))


def error_norms_data(initial_condition, flux): 
 """!
 @brief Compute and store various norms of the error vector for all grid resolutions for given initial condition and flux.
 
 @param initial_condition One of the valid initial conditions. @see initial_conditions

 @param flux One of the valid fluxes. @see fluxes
 
 @return Array of dictionaries indexed by grid refinement exponent containing "norm : error" key-value pairs
 """
 
 # Path to the hdf5 file containing results of computations
 computations_database = h5py.File(database_path, "r")
  
# Array of dictionaries containing errors in various norms
 _error_norms = {}

# Compute and store various norms of the error vector for all grid resolutions (i.e. from k = 6 to k = 16)
 for k in range(6, 17):
# Paths to the relevant datasets   
  dataset_initial_path = initial_condition + "/k = " + str(k) + " initial_data"
  dataset_path = initial_condition + "/" + flux + "/k = " + str(k)
   
# Store the computational data as arrays to enable vector subtraction   
  initial_data = numpy.array(computations_database[dataset_initial_path])
  computed_data = numpy.array(computations_database[dataset_path])

# Error vector  
  error = (initial_data - computed_data)
# Store various norms of the error vector  
  _error_norms[k] = {
  "sup_norm" : grid_norm(error, numpy.inf), 
  "one_norm" : grid_norm(error, 1), 
  "two_norm" : grid_norm(error, 2)
  }
  
 computations_database.close()
 return _error_norms
 
def convergence_table_row_data(k, _error_norms):
 """!
 @brief Compute and store various norms of the error vector for all grid resolutions for given initial condition and flux.

 @param k Refinement exponent in the expression for grid stepsize: h = 2^(-k)
 
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
 data["k"] = k
 data["M"] = 2**k 
 data["h"] = 2**(-k) 
 data["sup_norm"] = _error_norms[k]["sup_norm"] 
 data["sup_norm_rate"] = convergence_rate(k, "sup_norm") 
 data["one_norm"] = _error_norms[k]["one_norm"]
 data["one_norm_rate"] = convergence_rate(k, "one_norm")
 data["two_norm"] = _error_norms[k]["two_norm"]
 data["two_norm_rate"] = convergence_rate(k, "two_norm")
 return data

def computation_attributes(initial_condition, flux):  
 """!
 @brief Retreive attributes pertaining to the computation.

 @param initial_condition One of the valid initial conditions. @see initial_conditions

 @param flux One of the valid fluxes. @see fluxes

 @return Dictionary of attributes pertaining to the computation.
  """ 
 
 group_path = initial_condition + "/" + flux 
 
 computations_database = h5py.File(database_path, "r")
 attributes = {}
 attributes["a"] = computations_database[group_path].attrs["a"]
 attributes["T"] = computations_database[group_path].attrs["T"]
 attributes["CFL"] = computations_database[group_path].attrs["CFL"]
 computations_database.close()
 return attributes

def create_convergence_table(initial_condition, flux):
 """! 

 @brief Form and write the table containing information about the convergence of the numerical method for the given initial condition.

 @param initial_condition One of the valid initial conditions. @see initial_conditions

 @param flux One of the valid fluxes. @see fluxes

 """

 group_path = initial_condition + "/" + flux 
 
 error_norms = error_norms_data(*group_path.split("/"))

# Return a formatted string containing convergence information for the given grid resolution 
 def convergence_row(k): 
   return "{k:^3} {M:^10} {h:^15.10f} {sup_norm:^15.10f} {sup_norm_rate:^5.2f} {one_norm:^15.10f} {one_norm_rate:^5.2f} {two_norm:^15.10f} {two_norm_rate:^5.2f}".format(**convergence_table_row_data(k, error_norms))

# Formatted string containing the text of the header of the table
 header = "{k:^3} {M:^10} {h:^15} {sup_norm:^15} {sup_norm_rate:^5} {one_norm:^15} {one_norm_rate:^5} {two_norm:^15} {two_norm_rate:^5}".format(**convergence_table_header_data)

# Formatted string containing the attributes of the computation
 closing = "Speed: {a_text:>12} {a} \nCFL Number: sigma = {CFL} \nOutput Time: {T_text:>6} {T}".format(
 a_text = "a =",
 T_text = "T =",
 **computation_attributes(*group_path.split("/"))
 )

# Write the table to a text file  
 output_text_table_file = open(html_path + "convergence_tables/" + group_path.replace("/", "_") + ".txt", "w") 
 output_text_table_file.write(header + "\n\n")
 for k in range(7, 17):
  output_text_table_file.write(convergence_row(k) + "\n")
 output_text_table_file.write("\n" + closing)

def output_html_content(initial_condition, flux):
 """! 

 @brief Form and write an HTML document containing:
 <ul>
  <li>
    plot of the computed solutions and the initial data
  </li>
  <li>
    convergence table
  </li>
  <li>
    analysis of the results
  </li>
 </ul>

 @param initial_condition One of the valid initial conditions. @see initial_conditions

 @param flux One of the valid fluxes. @see fluxes

 """  
 
 group_path = initial_condition + "/" + flux 
 
 error_norms = error_norms_data(*group_path.split("/"))

# Convergence table header represented as a list in order to iterate over its elements in the order the keys appear below.
# It is necessary because iterating over dictionary keys is done in random order, which will write the elements of 
# the table header out of order.
 dictionary_keys_ordered = [
  "k", 
  "M", 
  "h", 
  "sup_norm", 
  "sup_norm_rate", 
  "one_norm", 
  "one_norm_rate", 
  "two_norm", 
  "two_norm_rate"
 ]

# The format of the data to be included in the convergence table, represented as a dictionary in which the keys are 
# the names of the columns in the table
 html_table_format = {
  "k" : "{k}", 
  "M" : "{M}",
  "h" : "{h:.10f}",
  "sup_norm" : "{sup_norm:.10f}",
  "sup_norm_rate" : "{sup_norm_rate:.2f}",
  "one_norm" : "{one_norm:.10f}",
  "one_norm_rate" : "{one_norm_rate:.2f}",
  "two_norm" : "{two_norm:.10f}",
  "two_norm_rate" : "{two_norm_rate:.2f}"
 }

# Attributes of the various HTML elements:
# Convergence table belongs to the class "colored," which is used to control its appearance through CSS
 table_attributes = {"class" : "colored"}
# Image of the plot has attributes "src," which is the path to the image file, and it belong to "center" class,
# which means the image will be centered through CSS
 img_attributes = {"src" : "images/{0}.png".format(group_path.replace("/", "_")), "class" : "center"}
# Table caption is centered (the id "caption" is used in CSS file to control its appearance), and spans all the columns of the table.
 table_caption_attributes = {"colspan" : "9", "id" : "caption"}

# Write the name of the initial condition as a header (it will only be written preceding the first method "Upwind")
 h2 = etree.Element("h2")
# Link (anchor) to the element in the document: needed for acces drom the navigation menu  
 etree.SubElement(h2, "a", id = initial_condition.lower()).text = initial_condition

# Write the name of the numerical method as a header 
 h3 = etree.Element("h3")
# Link (anchor) to the element in the document: needed for acces drom the navigation menu
 etree.SubElement(h3, "a", id = group_path.replace("/", "_").lower()).text = flux
 
# Include the plot of the computations 
 figure = etree.Element("figure")
 
# Link to the plot in higher resolution: accessed by clicking on the image from the browser 
 a = etree.SubElement(figure, "a", href = "pdf/plots/{0}.pdf".format(group_path.replace("/", "_")))
 img = etree.SubElement(a, "img", **img_attributes)
# Caption of the plot 
 figcaption = etree.SubElement(figure, "figcaption").text = group_path.replace("/", ", ") + ": convergence visualization"

# Element corresponding to the convergence table in the document
 table = etree.Element("table")
# Element corresponding to the row of the convergence table containing its caption
 tr = etree.SubElement(table, "tr")
# Element corresponding to the header (caption) of the table
 th = etree.SubElement(tr, "th", **table_caption_attributes)
# Link to the text file containing convergence table 
 etree.SubElement(th, "a", href = "convergence_tables/{0}.txt".format(group_path.replace("/", "_"))).text = "Convergence Table" 

# Element corresponding to the header row of the convergence table
 tr = etree.SubElement(table, "tr")
 for key in dictionary_keys_ordered:
# Element of the header (e.g "Sup-Norm")
  etree.SubElement(tr, "th").text = convergence_table_header_data[key]

# Write the convergence data into the table
 for k in range(7,17):
# Element corresponding to the data row of the convergence table   
  tr = etree.SubElement(table, "tr") 
  for key in dictionary_keys_ordered:
# Color columns with headers "k" and "Rate" and write convergence data corresponding to the current row and column of the table
   if key == "k" or "rate" in key:
    etree.SubElement(tr, "td", **table_attributes).text = html_table_format[key].format(**convergence_table_row_data(k, error_norms))
# Write convergence data corresponding to the current row and column of the table    
   else:
    etree.SubElement(tr, "td").text = html_table_format[key].format(**convergence_table_row_data(k, error_norms))

# Display the text of the analysis of the computational results   
# content_id identifies the results of the computation in the HTML document 
 content_id = initial_condition + "_" + flux
# content_text_id identifies the text of the analysis of the computation  
 content_text_id = content_id + "_text"
# Block corresponding to the text of the analysis
 div_text = etree.Element("div", id=content_text_id)
# Script to include the text of the analysis 
 etree.SubElement(div_text, "script").text = "include_text(\"" + content_id + "\")"

# Write the resulting HMTL document into a file
 output_html_table_file = open(html_path + "html/computations/" + group_path.replace("/","_") + ".html", "w")
# Write the name of the initial condition as a header only preceding the first method ("Upwind") 
 if flux == "Upwind":
  output_html_table_file.write(etree.tostring(h2, pretty_print=True).decode("utf-8") + "\n")
 output_html_table_file.write(etree.tostring(h3, pretty_print=True).decode("utf-8") + "\n")
 output_html_table_file.write(etree.tostring(figure, pretty_print=True).decode("utf-8") + "\n")
 output_html_table_file.write(etree.tostring(table, pretty_print=True).decode("utf-8"))
 output_html_table_file.write(etree.tostring(div_text, pretty_print=True).decode("utf-8"))
 
def include_content_html():
 """! 

 @brief Create an HTML file responsible for the dynamical display of content (plots, convergence tables) for all computations in the website main page.
 
 Its purpose is to avoid inserting content into the HTML document manually, since it is highly error prone due to large amount og text to be manipulated.

 """

# <div id="content"> 
 root_div = etree.Element("div", id="content")
 
 for initial_condition in initial_conditions:
  for flux in fluxes:
   # content_id identifies the results of a particular computation in the HTML document 
   content_id = initial_condition + "_" + flux
   # <div id="content_id">
   div = etree.SubElement(root_div, "div", id=content_id)
   # JQuery function to include content dynamically
   # <script> = include_content(content_id)</script>
   etree.SubElement(div, "script").text = "include_content(\"" + content_id + "\")"
   #</div>   
# </div>

# Write the generated HTML document to a file
 output_html_file = open(html_path + "html/computations/include_content.html", "w") 
 output_html_file.write(etree.tostring(root_div, pretty_print=True).decode("utf-8"))