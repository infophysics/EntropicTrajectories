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
    
  }

}
