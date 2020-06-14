#include <iostream>
#include <string>
#include "vector.h"
#include "matrix.h"
#include "matrix.cpp"
#include "vector.cpp"
#include <lapacke.h>

int main()
{
  ET::Matrix<double> mat1("A",4,4,1.0);
  ET::Matrix<double> mat2("B",4,4,3.0);
  ET::Matrix<double> mat3 = mat1 + mat2;
  mat3.set_name("C");
  mat1.print();
  mat3.print();
  /* Locals */
    double A[5][3] = {1,1,1,2,3,4,3,5,2,4,2,5,5,4,3};
    double b[5][2] = {-10,-3,12,14,14,12,16,16,18,16};
    lapack_int info,m,n,lda,ldb,nrhs;

    /* Initialization */
    m = 5;
    n = 3;
    nrhs = 2;
    lda = 3;
    ldb = 2;



    /* Executable statements */
    printf( "LAPACKE_dgels (row-major, high-level) Example Program Results\n" );
    /* Solve least squares problem*/
    info = LAPACKE_dgels(LAPACK_ROW_MAJOR,'N',m,n,nrhs,*A,lda,*b,ldb);

    /* Print Solution */
    printf( "\n" );
  return 0;
}
