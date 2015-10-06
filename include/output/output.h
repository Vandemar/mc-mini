#pragma once

#include <fstream>

#include "geometry/geometry.h"
#include "problem/problem.h"
#include "parser/parser.h"

class OutputStructure {
  public:
    OutputStructure (ParamParser&       pp,
                     GeometryStructure& gs,
                     ProblemStructure&  ps);

    ~OutputStructure();

    void writeDefaultFile();
    void writeDefaultFile (const int timestep);

    void writeHDF5File();
    void writeHDF5File (int timestep);
  
  private:
    ParamParser&       parser;
    GeometryStructure& geometry;
    ProblemStructure&  problem;

    //The number of cells in the x-direction.
    int M; 
    //The number of cells in the y-direction.
    int N;

    double dx;

    string outputFormat;
    string outputPath;
    string outputFilename;

    std::ofstream problemXdmfFile;
};
