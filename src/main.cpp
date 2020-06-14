#include <iostream>
#include <string>
#include <vector>
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

  ET::Matrix<double> mat4("D",2,{1.0,2.0,3.0,4.0});
  ET::Matrix<double> mat5("E",2,{-1.0,2.0,-3.0,4.0});
  mat4.print();
  mat5.print();
  ET::Matrix<double> mat6 = mat4 * mat5;
  mat6.set_name("F");
  mat6.print();

  ET::Matrix<double> mat7("G",2,3,{1,2,3,4,5,6});
  mat7.print();

  
  return 0;
}
