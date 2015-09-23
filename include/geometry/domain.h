#pragma once

#include "parser/parser.h"
#include "geometry/geometry.h"

class FieldSet{
  public:
    FieldSet();
    
    double uVelocity;
    double vVelocity; 
};

class Domain{
  public:
    Domain(GeometryStructure &geometry);
    ~Domain();
    FieldSet *domain;
  
  private:
    // Dimension of domain is M x N.
    int M; 
    int N;
};
