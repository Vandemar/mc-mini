#include <iostream>

#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

MatrixXd makeDivX (const unsigned int);
MatrixXd makeDivY (const unsigned int);
MatrixXd makeGradX (const unsigned int);
MatrixXd makeGradY (const unsigned int);
MatrixXd makeLaplacian (const unsigned int);

int main () {
  const unsigned int N = 4;

  const double viscosity = 1.0;

  double * rhsData = new double[(3 * N + 2) * N];
  double * uData = rhsData;
  double * vData = rhsData + (N * (N + 1));
  double * pData = rhsData + (N * (2 * N + 2));
  
  Map<VectorXd> rhsVect (rhsData, (3 * N + 2) * N);

  Map<VectorXd> uVect (uData, N * (N + 1));
  Map<Matrix<double, Dynamic, Dynamic, RowMajor> > uMatrix (uData, N, N + 1);

  Map<VectorXd> vVect (vData, N * (N + 1));
  Map<Matrix<double, Dynamic, Dynamic, RowMajor> > vMatrix (vData, N + 1, N);

  Map<VectorXd> pVect (pData,     N * N);
  Map<Matrix<double, Dynamic, Dynamic, RowMajor> > pMatrix (pData, N, N);

  MatrixXd A = MatrixXd::Zero (N * (N - 1) + N * (N - 1) + N * N, N * (N + 1) + N * (N + 1) + N * N);
  
  A.block (0, 0, N * (N - 1), N * (N + 1))                     = -viscosity * makeLaplacian (N); /* MatrixXd::Constant (N * (N - 1), N * (N + 1), 1); */
  A.block (N * (N - 1), N * (N + 1), N * (N - 1), N * (N + 1)) = -viscosity * makeLaplacian (N); /* MatrixXd::Constant (N * (N - 1), N * (N + 1), 1); */
  A.block (0, 2 * N * (N + 1), N * (N - 1), N * N)             = makeGradX (N); /* MatrixXd::Constant (N * (N - 1), N * N, 2); */
  A.block (N * (N - 1), 2 * N * (N + 1), N * (N - 1), N * N)   = makeGradY (N); /* MatrixXd::Constant (N * (N - 1), N * N, 3); */
  A.block (2 * N * (N - 1), 0, N * N, N * (N + 1))             = makeDivX (N); /* MatrixXd::Constant (N * N, N * (N + 1), 4); */
  A.block (2 * N * (N - 1), N * (N + 1), N * N, N * (N + 1))   = makeDivY (N); /* MatrixXd::Constant (N * N, N * (N + 1), 5); */

  uVect = VectorXd::LinSpaced (N * (N + 1), 0, N * (N + 1) - 1);
  vVect = VectorXd::LinSpaced (N * (N + 1), 0, N * (N + 1) - 1);
  pVect = VectorXd::LinSpaced (N * N,       0, N * N - 1);

  cout << A << endl;

  //  cout << pMatrix.colwise().reverse() << endl;
  
  return 0;
}

MatrixXd makeLaplacian (const unsigned int N) {
  MatrixXd laplacian = MatrixXd::Zero (N * (N - 1), N * (N + 1));
  MatrixXd laplacianBlock = MatrixXd::Zero (N - 1, N + 1);

  laplacianBlock.diagonal (1) = VectorXd::Constant (N - 1,  4);
  laplacianBlock.diagonal () = laplacianBlock.diagonal (2) = VectorXd::Constant (N - 1, -1); 
  laplacianBlock (0, 0) = laplacianBlock (N - 2, N) = 0;

  for (int i = 0; i < N; ++i)
    laplacian.block (i * (N - 1), i * (N + 1), (N - 1), (N + 1)) = laplacianBlock;

  laplacianBlock = MatrixXd::Zero (N - 1, N + 1);
  laplacianBlock.diagonal (1) = VectorXd::Constant (N - 1, -1);

  for (int i = 0; i < N - 1; ++i) 
    laplacian.block (i * (N - 1), (i + 1) * (N + 1), N - 1, N + 1) = laplacian.block ((i + 1) * (N - 1), i * (N + 1), N - 1, N + 1) = laplacianBlock;

  laplacian (0, 1) = laplacian (N * (N - 1) - 1, N * (N + 1) - 2) = 5;

  return laplacian;
}

MatrixXd makeDivX (const unsigned int N) {
  MatrixXd divX = MatrixXd::Zero (N * N, (N + 1) * N);
  MatrixXd divBlock = MatrixXd::Zero (N , N + 1);

  divBlock.diagonal () = VectorXd::Constant  (N, -1);
  divBlock.diagonal (1) = VectorXd::Constant (N,  1);
  
  // Handle zero boundaries by setting x-boundaries to 0.
  divBlock (0, 0) = divBlock (N - 1, N) = 0;

  for (int i = 0; i < N; ++i)
    divX.block (i * N, i * (N + 1), N, N + 1) = divBlock;  
    
  return divX;
}

MatrixXd makeDivY (const unsigned int N) {
  MatrixXd divY = MatrixXd::Zero (N * N, N * (N + 1));
  MatrixXd divBlock = MatrixXd::Zero (N, N);

  divY.diagonal ()  = VectorXd::Constant (N * N, -1);
  divY.diagonal (N) = VectorXd::Constant (N * N,  1);

  // Handle zero boundaries by setting y-boundaries to 0.
  divY.block (0, 0, N, N) = divY.block (N * (N - 1), N * N, N, N) = divBlock;

  return divY;
}

MatrixXd makeGradX (const unsigned int N) {
  MatrixXd gradX = MatrixXd::Zero (N * (N - 1), N * N);
  MatrixXd gradBlock = MatrixXd::Zero (N - 1, N);

  gradBlock.diagonal () = VectorXd::Constant  (N - 1, -1);
  gradBlock.diagonal (1) = VectorXd::Constant (N - 1,  1);

  for (int i = 0; i < N; ++i)
    gradX.block (i * (N - 1), i * N, N - 1, N) = gradBlock;
  
  return gradX;
}

MatrixXd makeGradY (const unsigned int N) {
  MatrixXd gradY = MatrixXd::Zero (N * (N - 1), N * N);
  
  gradY.diagonal ()  = VectorXd::Constant (N * (N - 1), -1);
  gradY.diagonal (N) = VectorXd::Constant (N * (N - 1),  1);
  
  return gradY;
}
