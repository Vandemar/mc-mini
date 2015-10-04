import lib_output_processing as output
import subprocess

output.test()

## Create HDF5 databse if it does not exist
#output.create_output_database()

## Empty task list
#task = ["", ""]

## Run computations for the selected initial conditions and grid resolutions
#for initial_condition in output.initial_conditions:
 #for flux in output.fluxes:
  #task[0] = initial_condition
  #task[1] = flux
  #for k in range(6,17):
   #subprocess.call(["./compute_task", task[0], task[1], str(k)])

## Process the results for the selected initial conditions and grid resolutions
#for initial_condition in output.initial_conditions:
 #for flux in output.fluxes:
  #task[0] = initial_condition
  #task[1] = flux
  
  #output.check_conservation(*task)
  #output.create_plots(*task)
  #output.create_convergence_table(*task)
  #output.output_html_content(*task)
  #output.include_content_html()