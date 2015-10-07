#include <iostream>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <Eigen/Sparse>

#include "matrixForms/sparseForms.h"
#include "geometry/dataWindow.h"
#include "geometry/geometry.h"
#include "problem/problem.h"
#include "parser/parser.h"
#include "debug.h"

/*
 *
 * Data initialization routines
 *
 */ 

void ProblemStructure::initializeProblem() {
  initializeTimestep();
  initializeTemperature();
  initializeTemperatureBoundary();
  initializeVelocityBoundary();
  initializeViscosity();
}

void ProblemStructure::initializeTimestep() {
  deltaT = cfl * h / diffusivity;
  int nTimestep = (endTime - time) / deltaT;
  if (abs (nTimestep * deltaT + time - endTime) > 1E-06) 
    deltaT = (endTime - time) / ++nTimestep;
  #ifdef DEBUG
    std::cout << "<Timestep initialized to " << deltaT << ">" << std::endl;
  #endif
}

void ProblemStructure::initializeTemperature() {
  DataWindow<double> temperatureWindow (geometry.getTemperatureData(), M, N);

  double referenceTemperature;
  double temperatureScale;

  if (parser.push ("problemParams")) {
    if (parser.push ("initialTemperatureParams")) {
      parser.queryParamDouble ("referenceTemperature", referenceTemperature, 273.15);
      parser.queryParamDouble ("temperatureScale",     temperatureScale,     100.0);

      parser.pop();
    } 
    parser.pop();
  }
  
  if (temperatureModel == "constant") {
    
    for (int j = 0; j < N; ++j)
      for (int i = 0; i < M; ++i)
        temperatureWindow (i, j) = referenceTemperature;
   
  } else if (temperatureModel == "sineWave") {
    int xModes;
    int yModes;

    if (parser.push ("problemParams")) {
      if (parser.tryPush ("initialTemperatureParams")) {
        parser.queryParamInt ("xModes", xModes, 2);
        parser.queryParamInt ("yModes", yModes, 2);

        parser.pop();
      } 

      parser.pop();
    }

    for (int j = 0; j < N; ++j)
      for (int i = 0; i < M; ++i)
        temperatureWindow (i, j) = referenceTemperature +
                                   sin ((i + 0.5) * h * xModes * M_PI / xExtent) * 
                                   sin ((j + 0.5) * h * yModes * M_PI / yExtent) * 
                                   temperatureScale;

  } else if (temperatureModel == "squareWave") {
    //TODO: Swap the indicies. -HL
    for (int j = 0; j < N; ++j) 
      for (int i = 0; i < M; ++i) {
        if ((xExtent * 0.25 <= i*h && i*h <= xExtent * 0.75) && (yExtent * 0.25 <= j*h && j*h <= yExtent * 0.75))
          temperatureWindow (i, j) = referenceTemperature + temperatureScale;
        else
          temperatureWindow (i, j) = referenceTemperature;
      }
  } else if (temperatureModel == "fallingSquare") {

   /* In this implementation, the temperature functions merely as a compositional field with a value
      of 1 inside the square and 0 outside the square.
      The ratio of the area of the falling square to the grid area of the rectangular domain is 1 to 25.
      The square 100 km x 100 km = 100000 m X 100000 m is placed in the region
      [200 km, 300 km] X [50 km, 150 km] = [200000 m, 300000 m] X [50000 m, 150000 m]; i.e.,
//ToDo Fix 0.4 <= dx <= 0.6 and 0.1 <= dy <= 0.3
   */

    for (int j = 0; j < N; ++j) {
      for (int i = 0; i < M; ++i) {
        if ((xExtent * 0.4 <= i*h && i*h <= xExtent * 0.6) && (yExtent * 0.7 >= j*h  && j*h <= yExtent * 0.9))
          temperatureWindow (i, j) = referenceTemperature + temperatureScale;
        else
          temperatureWindow (i, j) = referenceTemperature;
      }
    }  
  } else if (temperatureModel == "circle") {
     double center_x, center_y, radius;

     if (parser.push ("problemParams")) {
       if (parser.tryPush ("initialTemperatureParams")) {
         parser.getParamDouble ("radius", radius);
         parser.getParamDouble ("xCenter", center_x);
         parser.getParamDouble ("yCenter", center_y);
         parser.pop();
       }
       parser.pop();
     }

     for (int i = 0; i < M; ++i)
       for (int j= 0; j < N; ++j) {
         if ( std::sqrt(std::pow((i*h+h/2)-(xExtent * center_x),2.0) + std::pow((j*h+h/2)-(yExtent * center_y),2.0))  < 2*h )
           temperatureWindow (i, j) = referenceTemperature + temperatureScale;
         else
           temperatureWindow (i, j) = referenceTemperature; 
       }
  } else {
    cerr << "<Unexpected temperature model: \"" << boundaryModel << "\" : Shutting down now>" << endl;
    exit(-1);
  }

  #ifdef DEBUG
    cout << "<Initialized temperature model as: \"" << temperatureModel << "\">" << endl;
    cout << "<Temperature Data>" << endl;
    cout << temperatureWindow.displayMatrix() << endl << endl;
  #endif
}

void ProblemStructure::initializeTemperatureBoundary() {
  //TODO: Is the temperature boundary on the top and bottom or right and left of the grid?
  DataWindow<double> temperatureBoundaryWindow (geometry.getTemperatureBoundaryData(), M, 2);

  double upperTemperature;
  double lowerTemperature;

  if (parser.push ("problemParams")) {
    if (parser.push ("temperatureBoundaryParams")) {
      parser.getParamDouble ("upperBoundaryTemperature", upperTemperature);
      parser.getParamDouble ("lowerBoundaryTemperature", lowerTemperature);

      parser.pop();
    }
    
    parser.pop();
  }

  for (int i = 0; i < M; ++i) {
    temperatureBoundaryWindow (i, 0) = lowerTemperature;
    temperatureBoundaryWindow (i, 1) = upperTemperature;
  }
}

void ProblemStructure::initializeVelocityBoundary() {
  DataWindow<double> uVelocityBoundaryWindow (geometry.getUVelocityBoundaryData(), 2, N);
  DataWindow<double> vVelocityBoundaryWindow (geometry.getVVelocityBoundaryData(), M, 2);


  if (boundaryModel == "tauBenchmark") {
    //TODO: Fix the indexing scheme
    for (int j = 0; j < N; ++j)
      for (int i = 0; i < 2; ++i)
        uVelocityBoundaryWindow (i, j) = cos (i * M * h) * sin ((j + 0.5) * h);
    for (int j = 0; j < 2; ++j)
      for (int i = 0; i < M; ++i)
        vVelocityBoundaryWindow (i, j) = -sin ((i + 0.5) * h) * cos (j * N * h);
  } else if (boundaryModel == "solCXBenchmark" ||
             boundaryModel == "solKZBenchmark" ||
             boundaryModel == "noFlux") {
    for (int j = 0; j < N; ++j)
      for (int i = 0; i < 2; ++i)
        uVelocityBoundaryWindow (i, j) = 0;
    for (int j = 0; j < 2; ++j)
      for (int i = 0; i < M; ++i)
        vVelocityBoundaryWindow (i, j) = 0;
  } else {
    cerr << "<Unexpected boundary model: \"" << boundaryModel << "\" : Shutting down now>" << endl;
    exit(-1);
  }

  
  #ifdef DEBUG
    cout << "<Initialized boundary model as: \"" << boundaryModel << "\">" << endl;
    cout << "<U Velocity Boundary Data>" << endl;
    cout << uVelocityBoundaryWindow.displayMatrix() << endl;
    cout << "<V Velocity Boundary Data>" << endl;
    cout << vVelocityBoundaryWindow.displayMatrix() << endl << endl;
  #endif
}

void ProblemStructure::initializeViscosity() {
  DataWindow<double> viscosityWindow (geometry.getViscosityData(), M + 1, N + 1);

  double viscosity;

  if (viscosityModel == "constant") {
    if (parser.push ("problemParams")) {
      if (parser.tryPush ("initialViscosity")) {
        parser.queryParamDouble ("viscosityScale", viscosity, 1.0E21);
      
        parser.pop();
      } else {
        viscosity = 1.0E21;
      }
      parser.pop();
    }

    for (int j = 0; j < (N + 1); ++j)
      for (int i = 0; i < (M + 1); ++i)
        viscosityWindow (i, j) = viscosity;
  } else if (viscosityModel == "tauBenchmark") {
    viscosity = 1.0;
  } else if (viscosityModel == "solCXBenchmark") {
    for (int j = 0; j < (N + 1); ++j) 
      for (int i = 0; i < (M + 1); ++i)
        viscosityWindow (i, j) = (i <= M / 2) ? 1.0 : 1.0E06;
  } else if (viscosityModel == "solKZBenchmark") {
    for (int j = 0; j < (N + 1); ++j)
      for (int i = 0; i < (M + 1); ++i)
        viscosityWindow (i, j) = 1.0 + i * h * 1.0E06;
  } else {
    cerr << "Unexpected viscosity model: \"" << viscosityModel << "\" : Shutting down now!" << endl;
    exit(-1);
  }

  #ifdef DEBUG
    cout << "<Viscosity model initialized as: \"" << viscosityModel << "\">" << endl;
    cout << "<Viscosity Data>" << endl;
    cout << viscosityWindow.displayMatrix() << endl << endl;
  #endif
}
