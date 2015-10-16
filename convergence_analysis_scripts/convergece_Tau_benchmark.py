from math import log, sqrt, exp, cos, sin, pi
import h5py
import itertools
import numpy


#Richardson:
#for i in range(0,[N_x/2])
#for j in range(0,[N_y/2])
#data_on_restricted_grid[i][j] = 0.25*(data[2*i][2*j]+data[2*i][2*j+1]+data[2*i+1][2*j]+data[2*i+1][2*j+1])

#Interpolation
  #for (int i = 0; i < M; ++i)
    #for (int j = 0; j < N; ++j)
      #if (j == 0) {
        #interpolatedUVelocityWindow (j, i) = (uVelocityBoundaryWindow (0, i) +
                                              #uVelocityWindow         (j, i)) / 2;
      #} else if (j == (N - 1)) {
        #interpolatedUVelocityWindow (j, i) = (uVelocityWindow         (j - 1, i) +
                                              #uVelocityBoundaryWindow (1, i)) / 2;
      #} else {
        #interpolatedUVelocityWindow (j, i) = (uVelocityWindow (j - 1, i) +
                                              #uVelocityWindow (j,     i)) / 2;
      #}

