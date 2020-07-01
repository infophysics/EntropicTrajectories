#include <iostream>
#include <string>
#include <vector>
#include "vector.h"
#include "matrix.h"
#include "ugrid.h"
#include "utils.h"
#include <lapacke.h>
#include <cblas.h>

int main()
{
  ET::Matrix<double> mat1("A",4,4,1.0);
  ET::UGrid<double> g("2D",2);

  std::cout << mat1.summary();

	ET::Matrix<double> id  = ET::identity_d(10);

	std::cout << id.summary();

	unsigned int N = 10000;
	unsigned int dim = 25;
	double max_range = 10.0;
	std::vector<std::vector<double> > samples;
	samples.resize(N);
  for (size_t i = 0; i < N; i++)
	{
		samples[i].resize(dim);
		for (size_t d = 0; d < dim; d++)
			samples[i][d] = max_range * (rand() % 1000) / (1000.0);
	}
	g.setUGrid(samples);
	g.queryNeighbors(10);



  return 0;
}
