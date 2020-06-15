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

int main()
{
  ET::Matrix<double> mat1("A",4,4,1.0);
  ET::Grid<double> g("2D",2);

  std::cout << mat1;

	unsigned int N = 100;
	unsigned int dim = 10;
	double max_range = 10.0;
	std::vector<std::vector<double> > samples;
	samples.resize(N);
  for (size_t i = 0; i < N; i++)
	{
		samples[i].resize(dim);
		for (size_t d = 0; d < dim; d++)
			samples[i][d] = max_range * (rand() % 1000) / (1000.0);
	}
	std::cout << "done\n";
	g.set_grid(samples);
	g.find_neighbors(3);

	

  return 0;
}
