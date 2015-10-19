#pragma once

#include <iostream>
#include <cassert>

#include <Eigen/Dense>

// Column-major data window per Eigen.

template<typename T>
class DataWindow {
  public:
    DataWindow (T* basePtr,
		unsigned int nXCells,
		unsigned int nYCells) :
		_basePtr(basePtr),
		_nXCells(nXCells), 
    _nYCells(nYCells)
		{};

    T& operator() (unsigned int x_i, unsigned int y_i) 
    {
      // Ensure we haven't gone out-of-bounds on memory. There may be cases
      // where we actually want to do that, but we can remove the assertion
      // if that actually happens.
      assert(y_i * _nXCells + x_i < _nXCells * _nYCells);
      return _basePtr[y_i * _nXCells + x_i];
    }
    
    const std::string displayMatrix() 
      {
	std::cout << Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> >(_basePtr, _nYCells, _nXCells).colwise().reverse();
  return "";
      }

  private:
    T *const           _basePtr;
    // Number of cells in the x-direction
    const unsigned int _nXCells;
    // Number of cells in the y-direction
    const unsigned int _nYCells;
};

