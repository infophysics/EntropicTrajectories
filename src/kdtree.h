//  Wrapper for a kDTree implementation
#pragma once

#include "matrix.h"
#include "alglibmisc.h"

namespace ET
{
  template<typename T>
  class kDTree
  {
  public:
    kDTree();
    ~kDTree();
    kDTree(Matrix<T> m);
  private:
    unsigned int _dim;
    unsigned int _N;
    
  };
}
