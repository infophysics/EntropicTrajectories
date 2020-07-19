#include "kdtree.h"
#include <random>
#include <iostream>
#include <memory>

using namespace ET;

void testKDTreeCreate()
{
  using array = std::vector<std::vector<double>>;
  //  generate some random data
  const int range_from  = 0.0;
  const int range_to    = 1.0;
  std::random_device rand_dev;
  std::mt19937 generator(rand_dev());
  std::uniform_real_distribution<double> distr(range_from, range_to);

  size_t N = 10000;
  std::shared_ptr<array> data = std::make_shared<array>(N);
  for (size_t i = 0; i < N; i++) {
   (*data)[i][0] = distr(generator);
   (*data)[i][1] = distr(generator);
  }

  KDTree<double> kdt(data);

  //  Query k neighbors for each point in data
  size_t k = 5;
  kdt.queryNeighbors(k);
}


int main(int, char**)
{
  void testKDTreeCreate();
  return 0;
}
