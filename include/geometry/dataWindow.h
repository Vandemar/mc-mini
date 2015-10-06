#pragma once

#include <string>
#include <iostream>
#include <Eigen/Dense>

// Column-major data window per Eigen.

template<typename T>
class DataWindow {
  public:
    DataWindow (T* basePtr = nullptr,
		unsigned int nXCells = 0,
		unsigned int nYCells = 0
	        ) :
		_basePtr(basePtr),
		_nXCells(nXCells), 
    _nYCells(nYCells)
		{};
	
  //Memory must be ordered in row-major format because hdf5 supports only row-major format.
  //I have had no luck finding a transpose funtion in hdf5. -HL
    T& operator() (unsigned int x_i, unsigned int y_i) 
      {
	return _basePtr[y_i * _nXCells + x_i];
      }
    
    void displayMatrix() 
      {
	std::cout << Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> >(_basePtr, _nYCells, _nXCells).colwise().reverse();
      }

  private:
    T *const           _basePtr;
    // Number of cells in the x-direction
    const unsigned int _nXCells;
    // Number of cells in the y-direction
    const unsigned int _nYCells;
};