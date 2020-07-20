//  speed_test_dot_product.cpp
#include "vector.h"
#include <random>
#include <fstream>
#include <iterator>
#include <memory>
#include <vector>
#include <chrono>

using namespace ET;

void speedTestDotProduct()
{
  //  The maximum number of elements to take the dot product of
  const int START = 100000;
  const int END = 10000000;
  const int INTERVAL = 10000;
  //  Arrays to hold the time differences
  std::vector<double> time_diff1;
  std::vector<double> time_diff2;
  std::vector<double> time_diff3;
  //  Setup vectors of coefficients
  std::vector<double> vec1;
  std::vector<double> vec2;
  //  Setup random number generator
  const int range_from  = 0.0;
  const int range_to    = 1.0;
  std::random_device rand_dev;
  std::mt19937 generator(rand_dev());
  std::uniform_real_distribution<double> distr(range_from, range_to);
  //  Setup clocks
  auto start = std::chrono::steady_clock::now();
  auto end = std::chrono::steady_clock::now();
  auto start2 = std::chrono::steady_clock::now();
  auto end2 = std::chrono::steady_clock::now();
  auto start3 = std::chrono::steady_clock::now();
  auto end3 = std::chrono::steady_clock::now();

  double duration = 0;
  double duration2 = 0;
  double duration3 = 0;
  double dot_brute = 0;
  double ddot = 0;
  double dot_cblas = 0;
  //  Iterate from 0 to NUM_ELEMENTS
  for(auto i = START; i < END; i+=INTERVAL) {
    //  Construct the two vectors using random elements
    vec1.resize(i);
    vec2.resize(i);
    for (auto j = 0; j < i; j++) {
      vec1[j] = distr(generator);
      vec2[j] = distr(generator);
    }
    Vector<double> v1(vec1);
    Vector<double> v2(vec2);
    //  Perform the brute force dot product
    // std::cout << "brute\n";
    //   Timer timer;
    start = std::chrono::steady_clock::now();
    dot_brute = v1.dot(v2);

    end = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>
                            (end-start).count();
    //  Store the result1
    time_diff1.push_back(duration);
    //  Perform the secont dot product
    start2 = std::chrono::steady_clock::now();
    // std::cout << "ddot\n;
    // Timer timer;
    ddot = DDOT(v1,v2);
    end2 = std::chrono::steady_clock::now();
    duration2 = std::chrono::duration_cast<std::chrono::nanoseconds>
                            (end2-start2).count();
    //  Store the result
    time_diff2.push_back(duration2);
    //  Try using cblas directly
    start3 = std::chrono::steady_clock::now();
    dot_cblas = cblas_ddot(i,//  dimension of the t_vectors
                      vec1.data(),  //  pointer to the elements of v
                      1,           //  increment of the elements of v
                      vec2.data(),  //  pointer to the elements of u
                      1);          //  increment of the elements of u

    end3 = std::chrono::steady_clock::now();
    duration3 = std::chrono::duration_cast<std::chrono::nanoseconds>
                            (end3-start3).count();
    //  Store the result1
    time_diff3.push_back(duration3);
  }
  //  Save the results to file
  {
    std::ofstream fout("speed_test_results_brute.txt");
    fout.precision(10);
    std::copy(time_diff1.begin(), time_diff1.end(),
              std::ostream_iterator<double>(fout, "\n"));
  }
  {
    std::ofstream fout("speed_test_results_ddot.txt");
    fout.precision(10);
    std::copy(time_diff2.begin(), time_diff2.end(),
              std::ostream_iterator<double>(fout, "\n"));
  }
  {
    std::ofstream fout("speed_test_results_cblas_direct.txt");
    fout.precision(10);
    std::copy(time_diff3.begin(), time_diff3.end(),
              std::ostream_iterator<double>(fout, "\n"));
  }
}

int main()
{
  speedTestDotProduct();
}
