#include "geometry/domain.h"

Domain::Domain(GeometryStructure &geometry) {
  // Originally the domain was swapped such that M is the number of columns and N is the number of rows. 
  M = geometry.getN();
  N = geometry.getM();
  domain = new FieldSet[M*N];
}

Domain::~Domain() {
  delete[] domain;
}

FieldSet::FieldSet() {
  uVelocity = vVelocity = 0;
}
