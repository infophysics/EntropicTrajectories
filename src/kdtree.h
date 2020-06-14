//  Wrapper for a KDTree implementation
#pragma once

#include "alglibmisc.h"

namespace ET
{
  template<typename T>
  class KDTree
  {
  public:
    KDTree();
    ~KDTree();
  private:
    alglib::kdtree _kdt;

  };
}
