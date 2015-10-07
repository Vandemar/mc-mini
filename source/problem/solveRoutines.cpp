#include <iostream>
#include <cassert>
#include <cmath>

#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "debug.h"

#ifndef USE_DENSE
#include "matrixForms/sparseForms.h"
#else
#include "matrixForms/denseForms.h"
#endif
#include "geometry/dataWindow.h"
#include "geometry/geometry.h"
#include "problem/problem.h"
#include "parser/parser.h"

using namespace Eigen;
using namespace std;

/*
 *
 * Main loop routines
 *
 */ 

// Update the forcing terms
// T -> F
void ProblemStructure::updateForcingTerms() {
  DataWindow<double> uForcingWindow (geometry.getUForcingData(), M - 1, N);
  DataWindow<double> vForcingWindow (geometry.getVForcingData(), M, N - 1);

  #ifdef DEBUG
    cout << "<Calculating forcing model using \"" << forcingModel << "\">" << endl;
  #endif

  if (forcingModel == "tauBenchmark") {
    // Benchmark taken from Tau (1991; JCP Vol. 99)
    for (int j = 0; j < N; ++j)
      for (int i = 0; i < M-1; ++i)
        uForcingWindow (i, j) = 3 * cos ((i + 1) * h) * sin ((j + 0.5) * h);

    for (int j = 0; j < N-1; ++j)
      for (int i = 0; i < M; ++i)
        vForcingWindow (i, j) = -sin ((i + 0.5) * h) * cos ((j + 1) * h);

  } else if (forcingModel == "solCXBenchmark" ||
             forcingModel == "solKZBenchmark") {
    // solCX Benchmark taken from Kronbichler et al. (2011)
    for (int j = 0; j < N; ++j)
      for (int i = 0; i < M-1; ++i)
        uForcingWindow (i, j) = 0;

    for (int j = 0; j < N-1; ++j)
      for (int i = 0; i < M; ++i)
        vForcingWindow (i, j) = - sin((j + 0.5) * M_PI * h) * cos ((i + 1) * M_PI * h);

  } else if (forcingModel == "vorticalFlow") {
    for (int j = 0; j < (N - 1); j++) 
      for (int i = 0; i < M; ++i)
        uForcingWindow (i, j) = cos ((i + 1) * h) * sin ((j + 0.5) * h);

    for (int j = 0; j < N-1; ++j) 
      for (int i = 0; i < M; ++i)
        vForcingWindow (i, j) = -sin ((i + 0.5) * h) * cos ((j + 1) * h);

  } else if (forcingModel == "buoyancy") {
    DataWindow<double> temperatureWindow (geometry.getTemperatureData(), M, N);

    double referenceTemperature;
    double densityConstant;
    double thermalExpansion;

    if (parser.push ("problemParams")) {
      if (parser.push ("buoyancyModelParams")) {
        parser.queryParamDouble ("referenceTemperature", referenceTemperature, 273.15);
        parser.queryParamDouble ("densityConstant",      densityConstant,      100.0);
        parser.queryParamDouble ("thermalExpansion",     thermalExpansion,       1.0);
        parser.pop();
      }
      parser.pop();
    }

    for (int j = 0; j < N; ++j)
      for (int i = 0; i < M-1; ++i)
        uForcingWindow (i, j) = 0;

// TODO: EGP thinks the average should be ((temperatureWindow(j,i)+temperatureWindow(j+1,i))/2)  
    for (int j = 0; j < N-1; ++j) {
      for (int i = 0; i < M; ++i)
        vForcingWindow (i, j) =  -1 * densityConstant *
                                  (1 - thermalExpansion * 
                                   ((temperatureWindow (i, j) + 
                                     temperatureWindow (i, j + 1)) / 2 -
                                      referenceTemperature));
      }
  } else if (forcingModel == "fallingSquare") {
    DataWindow<double> temperatureWindow (geometry.getTemperatureData(), N, M);

    double referenceTemperature;
    double densityConstant;
    double thermalExpansion;

    if (parser.push ("problemParams")) {
      if (parser.push ("buoyancyModelParams")) {
        parser.queryParamDouble ("referenceTemperature", referenceTemperature, 0.0);
        parser.queryParamDouble ("densityConstant",      densityConstant,   3200.0);
        parser.queryParamDouble ("thermalExpansion",     thermalExpansion,  -100.0);
        parser.pop();
      }
      parser.pop();
    }

    for (int j = 0; j < N; ++j)
      for (int i = 0; i < M-1; ++i)
        uForcingWindow (i, j) = 0;

    for (int j = 0; j < N-1; ++j) {
      for (int i = 0; i < M; ++i)
        vForcingWindow (i, j) =  -1 * densityConstant *
                                  (1 - thermalExpansion * 
                                   ((temperatureWindow (i, j) + 
                                     temperatureWindow (i+1, j)) / 2 -
                                      referenceTemperature));
      }
 } else {
    cerr << "<Unexpected forcing model: \"" << forcingModel << "\" : Shutting down now>" << endl;
    exit(-1);
  }

  #ifdef DEBUG
    cout << "<U Forcing Data>" << endl;
    cout << uForcingWindow.displayMatrix() << endl;
    cout << "<V Forcing Data>" << endl;
    cout << vForcingWindow.displayMatrix() << endl << endl;
  #endif
}

// Solve the stokes equation
// F -> U X P
void ProblemStructure::solveStokes() {
  Map<VectorXd> stokesSolnVector (geometry.getStokesData(), N * (M - 1) + (N - 1) * M + N * M);
  
  #ifndef USE_DENSE
  static SparseMatrix<double> stokesMatrix   (3 * N * M - N - M, 3 * N * M - N - M);
  static SparseMatrix<double> forcingMatrix  (3 * N * M - N - M, 2 * N * M - N - M);
  static SparseMatrix<double> boundaryMatrix (3 * N * M - N - M, 2 * N + 2 * M);
  static SparseLU<SparseMatrix<double>, COLAMDOrdering<int> > solver;
  #else
  /* Don't use this unless you hate your computer. */
  static MatrixXd stokesMatrix   (3 * M * N - M - N, 3 * M * N - M - N);
  static MatrixXd forcingMatrix  (3 * M * N - M - N, 2 * M * N - M - N);
  static MatrixXd boundaryMatrix (3 * M * N - M - N, 2 * M + 2 * N);
  static PartialPivLU<MatrixXd> solver;
  #endif

  static bool initialized;

  double * viscosityData = geometry.getViscosityData();

  if (!(initialized) || !(viscosityModel=="constant")) {

  #ifndef USE_DENSE
    SparseForms::makeStokesMatrix   (stokesMatrix,   N, M, h, viscosityData);
    stokesMatrix.makeCompressed();
    SparseForms::makeForcingMatrix  (forcingMatrix,  N, M);
    forcingMatrix.makeCompressed();
    SparseForms::makeBoundaryMatrix (boundaryMatrix, N, M, h, viscosityData);
    boundaryMatrix.makeCompressed();

    solver.analyzePattern (stokesMatrix);
    solver.factorize (stokesMatrix);
  #else
    DenseForms::makeStokesMatrix   (stokesMatrix, M, N, h, viscosityData);
    DenseForms::makeForcingMatrix  (forcingMatrix, M, N);
    DenseForms::makeBoundaryMatrix (boundaryMatrix, M, N, h, viscosityData);

    solver.compute (stokesMatrix);
  #endif
    initialized = true;
  }

  stokesSolnVector = solver.solve (forcingMatrix  * Map<VectorXd>(geometry.getForcingData(), 2 * N * M - N - M) + 
                                   boundaryMatrix * Map<VectorXd>(geometry.getVelocityBoundaryData(), 2 * N + 2 * M));

  Map<VectorXd> pressureVector (geometry.getPressureData(), N * M);
  double pressureMean = pressureVector.sum() / (N * M);
  pressureVector -= VectorXd::Constant (N * M, pressureMean);

#ifdef DEBUG
  cout << "<Calculated Stokes Equation Solutions>" << endl;
  cout << "<U Velocity Data>" << endl;
  cout << DataWindow<double> (geometry.getUVelocityData(), N - 1, M).displayMatrix() << endl;
  cout << "<V Velocity Data>" << endl;
  cout << DataWindow<double> (geometry.getVVelocityData(), N, M - 1).displayMatrix() << endl;
  cout << "<Pressure Data>" << endl;
  cout << DataWindow<double> (geometry.getPressureData(), N, M).displayMatrix() << endl << endl;
#endif
}

// Solve the advection/diffusion equation
// U X T -> T
void ProblemStructure::solveAdvectionDiffusion() {
  #ifdef DEBUG
    cout << "<Using \"" << advectionMethod << "\" for advection>" << endl;
  #endif
  if (advectionMethod == "upwindMethod") {
    upwindMethod();
  } else if (advectionMethod == "frommMethod") {
    frommMethod();
  } else if (advectionMethod == "none") {
  } else {
    cerr << "<Unexpected advection method: \"" << advectionMethod << "\" : Shutting down now>" << endl;
    exit (-1);
  }
  
  #ifdef DEBUG
    cout << "<Using \"" << diffusionMethod << "\" for diffusion>" << endl;
  #endif
  if (diffusionMethod == "forwardEuler") {
    forwardEuler();
  } else if (diffusionMethod == "backwardEuler") {
    backwardEuler();
  } else if (diffusionMethod == "crankNicolson") {
    crankNicolson();
  } else if (diffusionMethod == "none") {
  } else {  
    cerr << "<Unexpected diffusion method: \"" << diffusionMethod << "\" : Shutting down now>" << endl;
    exit (-1);
  }

  #ifdef DEBUG
    cout << "<Finished Advection/Diffusion Step>" << endl << endl;
  #endif
}
