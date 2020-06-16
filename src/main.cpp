#include <iostream>
#include <string>
#include <vector>
#include "vector.h"
#include "matrix.h"
#include "grid.h"
#include "matrix.cpp"
#include "vector.cpp"
#include "grid.cpp"
#include <lapacke.h>
#include <cblas.h>

int main()
{
  ET::Matrix<double> mat1("A",4,4,1.0);
  ET::Grid<double> g("2D",2);

  std::cout << mat1;

	ET::Matrix<double> id  = ET::identity<double>(10);

	std::cout << id;

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
	g.set_grid(samples);
	g.find_neighbors(10);



  return 0;
}
