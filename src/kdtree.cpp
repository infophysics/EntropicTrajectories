#include "kdtree.h"

namespace ET
{

  template<typename T>
  kDTree<T>::kDTree()
  {

  }
  template<typename T>
  kDTree<T>::~kDTree()
  {

  }
  template<typename T>
  kDTree<T>::kDTree(Matrix<T> m) : _N(m.get_rows()), _dim(m.get_cols())
  {
    _data.setcontent(_N, _dim, & *m.get_mat().begin());
    alglib::kdtreebuild(_data ,_N, _dim, 0, 2, _kdt);
    alglib::real_1d_array x = "[3,2]";
    alglib::ae_int_t k = alglib::kdtreequeryknn(_kdt,x,1);
    std::cout << k << std::endl;
  }

}
