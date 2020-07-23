#include "scalarfield.h"
#include "ugrid.h"
#include "interpolator.h"
#include "localtaylor.h"

#include <random>
#include <iostream>
#include <memory>


using namespace ET;

//  Test the functionality of derivatives
//  using the local Taylor expansion interpolator
//  on a scalar field in one dimension.
void testScalarFieldDerivativesLTE1D()
{
  //  Generate a uniform random sample space
  using array = std::vector<std::vector<double>>;
  //  generate some random data
  const int range_from  = -3.0;
  const int range_to    =  3.0;
  std::random_device rand_dev;
  std::mt19937 generator(rand_dev());
  std::uniform_real_distribution<double> distr(range_from, range_to);

  size_t N = 10000;
  array data(N);

  for (size_t i = 0; i < N; i++) {
    double x = distr(generator);
    data[i][0] = x;
  }

  //  Create the unstructured grid
  UGrid<double> ugrid("x",data);
  //----------------------------------------------------------------------------
  //  From here we can test different one dimensional scalar field values using
  //  the same underlying grid.

  //  (1) f(x) = cos(x)
  std::vector<double> cos_x(N);
  for (size_t i = 0; i < N; i++) {
    cos_x[i] = cos(data[i][0]);
  }
  //  Now create the scalar field object
  ScalarField<double> f_1("cos(x)",ugrid,cos_x);
  //  Adjust the settings of the interpolator
  f_1.getInterpolator()->set_k(5);
  f_1.getInterpolator()->setSearchScheme(SearchScheme::NEAREST_NEIGHBORS);
  //  Compute the first five derivatives and save the results to a file
  for (auto i = 1; i < 6; i++) {
    array f_1x = f_1.derivative(i);
    std::string filename = "f_1" + std::to_string(i) + ".csv";
    writeToCSV(filename,f_1x,",");
  }

}


int main(int, char**)
{
  testScalarFieldDerivativesLTE1D();
}
