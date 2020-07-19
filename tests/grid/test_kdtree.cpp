#include "kdtree.h"
#include <random>
#include <iostream>
#include <memory>

using namespace ET;

void testKDTreeDefaultConstructor()
{
  KDTree<double> kdt();
}

void testKDTreeConstructors()
{
  using array = std::vector<std::vector<double>>;
  //  generate some random data
  const int range_from  = 0.0;
  const int range_to    = 1.0;
  std::random_device rand_dev;
  std::mt19937 generator(rand_dev());
  std::uniform_real_distribution<double> distr(range_from, range_to);

  size_t N = 10000;
  array data(N);

  for (size_t i = 0; i < N; i++) {
    double x = distr(generator);
    double y = distr(generator);
    data[i].resize(2);
    data[i][0] = x;
    data[i][1] = y;
  }
  KDTree<double> kdt(std::make_shared<array>(data));
}

int main(int, char**)
{
  testKDTreeDefaultConstructor();
  testKDTreeConstructors();
}
