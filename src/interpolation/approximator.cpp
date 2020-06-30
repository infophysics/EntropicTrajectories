//------------------------------------------------------------------------------
//  approximator.cpp
//  The Entropic Trajectories Framework
//  -----------------------------------
//  Copyright (C) [2020] by [N. Carrara, F. Costa, P. Pessoa]
//  [ncarrara@albany.edu,felipecosta.physics@gmail.com,
//    pedroh.pessoa100@gmail.com]
//
//  Permission to use, copy, modify, and/or distribute this software for any
//  purpose with or without fee is hereby granted.
//
//  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
//  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
//  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
//  SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
//  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
//  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR
//  IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
//------------------------------------------------------------------------------
#include "approximator.h"


namespace ET
{
  //----------------------------------------------------------------------------
  //  Enum maps
  //----------------------------------------------------------------------------
  //  map for a string to enum of approximator type
  std::map<std::string, ApproxType> ApproxTypeMap =
  {
    { "MLS", ApproxType::MLS },
    { "RBF", ApproxType::RBF }
  };
  //  map for a string to enum of LSDriver type
  std::map<std::string, LSDriver> LSDriverMap =
  {
    { "xGELS",  LSDriver::xGELS },
    { "xGELSY", LSDriver::xGELSY },
    { "xGELSD", LSDriver::xGELSD },
    { "xGELSS", LSDriver::xGELSS }
  };
  std::map<ApproxType, std::string> ApproxTypeNameMap =
  {
    { ApproxType::MLS, "Moving least squares" },
    { ApproxType::RBF, "Radial basis functions" }
  };
  //  map for a string to enum of LSDriver type
  std::map<LSDriver, std::string> LSDriverNameMap =
  {
    { LSDriver::xGELS,  "xGELS" },
    { LSDriver::xGELSY, "xGELSY" },
    { LSDriver::xGELSD, "xGELSD" },
    { LSDriver::xGELSS, "xGELSS" }
  };
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Constructors
  //----------------------------------------------------------------------------
  template<typename T>
  Approximator<T>::Approximator()
  {
    std::cout << "\nApproximator created at location " << this;
    _type = 0;
    _lsdriver = 0;
  }
  template<typename T>
  Approximator<T>::~Approximator()
  {
    std::cout << "\nApproximator at location " << this << " destroyed.";
  }
  template<typename T>
  Approximator<T>::Approximator(int type) : _type(type)
  {
    std::cout << "\nApproximator created at location " << this;
  }
  template<typename T>
  Approximator<T>::Approximator(std::string type)
  {
    std::cout << "\nApproximator created at location " << this;
    _type = ApproxTypeMap[type];
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Getters and Setters
  //----------------------------------------------------------------------------
  template<typename T>
  int Approximator<T>::getApproxType() const
  {
    return _type;
  }
  template<typename T>
  ApproxParams Approximator<T>::getApproxParams() const
  {
    return _params;
  }
  template<typename T>
  int Approximator<T>::getLSDriver() const
  {
    return _lsdriver;
  }
  template<typename T>
  int Approximator<T>::getFlag() const
  {
    return _flag;
  }
  template<typename T>
  std::string Approximator<T>::getInfo() const
  {
    return _info;
  }
  template<typename T>
  void Approximator<T>::setApproxType(std::string type)
  {
    _type = ApproxTypeMap[type];
  }
  template<typename T>
  void Approximator<T>::setApproxParams(ApproxParams params)
  {
    _params = params;
  }
  template<typename T>
  void Approximator<T>::setLSDriver(std::string type)
  {
    _lsdriver = LSDriverMap[type];
  }
  template<typename T>
  void Approximator<T>::set_k(uint64_t k)
  {
    _params.k = k;
  }
  template<typename T>
  void Approximator<T>::set_n(uint64_t n)
  {
    _params.n = n;
  }
  template<typename T>
  void Approximator<T>::setFlag(int flag)
  {
    _flag = flag;
  }
  template<typename T>
  void Approximator<T>::setInfo(std::string info)
  {
    _info = info;
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Methods for Scalar fields
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<T> Approximator<T>::scalarGradient(const std::shared_ptr<UGrid<T>> ugrid,
    const std::shared_ptr<ScalarField<T>> field, uint64_t index)
  {
    if (_type == 0)
      return scalarGradientMLS(ugrid, field, index);
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<T> Approximator<T>::scalarGradientMLS(const std::shared_ptr<UGrid<T>> ugrid,
    const std::shared_ptr<ScalarField<T>> field, uint64_t index)
  {
    std::vector<T> result;
    //  First, find the nearest neighbors associated to the point specified by
    //  index.
    ugrid->queryNeighbors(_params.k);
    std::vector<uint64_t> neighbors = ugrid->getNeighbors(index);
    //  Construct the matrix associated with the ugrid spacing
    Matrix<T> B = constructTaylorMatrix(ugrid, neighbors, index, _params.n);
    //  Construct the vector of field values associated to each point
    std::vector<T> field_neighbors(_params.k);
    for (uint32_t i = 0; i < _params.k; i++)
    {
      field_neighbors[i] = (*field)(neighbors[i]);
    }
    Vector<T> field_vals(field_neighbors);
    Vector<T> answer = DGELS(B,field_vals);
    return answer.getVec();
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //
  //----------------------------------------------------------------------------
  template<typename T>
  Matrix<T> Approximator<T>::constructTaylorMatrix(const std::shared_ptr<UGrid<T>> ugrid,
    const std::vector<uint64_t> neighbors, uint64_t index, uint64_t order)
  {
    std::cout << "\nneighbors location: " << &neighbors;
    std::cout << "\nGrid location: " << &ugrid;
    Monomial mono(ugrid->getDim(),order);
    std::vector<std::vector<double> > B;
    for (uint64_t i = 0; i < neighbors.size(); i++)
    {
      uint64_t id = neighbors[i];
      std::vector<double> temp = mono.taylorMonomialExpansion(ugrid->getPoint(index),
            ugrid->getPoint(id));
      B.push_back(temp);
    }
    Matrix<T> m("B",B);
    return m;
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Various functions
  //----------------------------------------------------------------------------
  template<typename T>
  std::string Approximator<T>::summary()
  {
    std::string s = "\nApproximator type: " + ApproxTypeNameMap[_type];
    if (_type == MLS)
    {
      s += "\nLeast squares driver type: " + LSDriverNameMap[_lsdriver];
    }
    s += "\nApproximator parameters - k = " + std::to_string(_params.k);
    s += "\n                          n = " + std::to_string(_params.n);
    return s;
  }
  //----------------------------------------------------------------------------
}
