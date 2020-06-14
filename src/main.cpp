#include <iostream>
#include <string>
#include <vector>
#include "vector.h"
#include "matrix.h"
#include "grid.h"
#include "kdtree.h"
#include "matrix.cpp"
#include "vector.cpp"
#include "grid.cpp"
#include <lapacke.h>

int main()
{
  ET::Matrix<double> mat1("A",4,4,1.0);
  ET::Grid<double> g("2D",2);

  std::cout << mat1;

  return 0;
}
