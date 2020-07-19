#include "vector.h"
#include <random>
#include <string>
#include <iostream>
#include <memory>
#include <vector>

using namespace ET;

void testVectorDefaultConstructor()
{
  Vector<double> v();
  Vector<double>* v_p = new Vector<double>();
  delete v_p;
}

void testVectorConstructors()
{
  size_t dim = 10;
  std::string name = "vector";
  std::vector<double> vec = {0,1,2,3,4,5,6};
  double init = 3.14159;

  Vector<double> v_dim(dim);
  Vector<double> v_name_dim(name,dim);
  Vector<double> v_vec(vec);
  Vector<double> v_name_vec(name,vec);
  Vector<double> v_dim_init(dim,init);
  Vector<double> v_name_dim_init(name,dim,init);
}

void testVectorGetters()
{
  Vector<double> v("name",{0,1,2,3,4});
  size_t dim = v.getDim();
  std::vector<double> vec = v.getVec();
  std::vector<double>* vec_p = v.accessVec();
  double* v_data = v.data();
  std::string name = v.getName();
  int flag = v.getFlag();
  std::string info = v.getInfo();
}

void testVectorSetters()
{
  Vector<double> v;
  v.setDim(10);
  v.setVec({0,1,2,3,4});
  v.setName("name");
  v.setFlag(1);
  v.setInfo("info");
}

void testVectorOperators()
{
  Vector<double> v("v",{1,2,3,4});
  Vector<double> u("u",{4,3,2,1});
  double s = 10.0;

  Vector<double> w = v;
  bool same = (v == u);
  bool not_same = (v != u);
  Vector<double> minus_v = -v;
  Vector<double> v_sum_u = v + u;
  v += u;
  Vector<double> v_sub_u = v - u;
  v -= u;
  double v_dot_u = v * u;
  double v_dot_u_dot = v.dot(u);
  Vector<double> v_plus_s = v + s;
  Vector<double> v_sub_s = v - s;
  Vector<double> v_times_s = v * s;
  Vector<double> v_div_s = v / s;
  v += s;
  v -= s;
  v *= s;
  v /= s;
  Vector<double> s_plus_v = s + v;
  Vector<double> s_sub_v = s - v;
  Vector<double> s_times_v = s * v;
  Vector<double> s_div_v = s / v;

  double val = v(1);
  v(1) = 5.0;
}

int main(int, char**)
{
  testVectorDefaultConstructor();
  testVectorConstructors();
  testVectorGetters();
  testVectorSetters();
  testVectorOperators();
}
