#pragma once

#include <iostream>
#include <cassert>

#include <Eigen/Dense>

// Column-major data window per Eigen.

template<typename T>
class DataWindow {
  public:
    DataWindow (T* basePtr = nullptr,
		unsigned int columns = 0,
		unsigned int rows = 0
	        ) :
		_basePtr(basePtr),
		_cols(columns),
		_rows(rows) 
		{};
		
    T& operator() (unsigned int col, unsigned int row) 
      {
	return _basePtr[row * _cols + col];
      }
    
    void displayMatrix() 
      {
	std::cout << Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> >(_basePtr, _rows, _cols).colwise().reverse();
      }

  private:
    T *const           _basePtr;
    const unsigned int _cols;
    const unsigned int _rows;
};

