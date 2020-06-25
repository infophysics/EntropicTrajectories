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
  //  map for a string to enum of derivative type
  std::map<std::string, ApproxType> ApproxTypeMap =
  {
    { "MLS", ApproxType::MLS },
    { "RBF", ApproxType::RBF }
  };

  template<typename T>
  Approximator<T>::Approximator()
  {

  }

  template<typename T>
  Approximator<T>::~Approximator()
  {

  }

  //  Setters
  template<typename T>
  void Approximator<T>::setDerivative(std::string type)
  {
    _type = ApproxTypeMap[type];
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
  std::vector<T> Approximator<T>::gradient(Grid<T>* grid,
    ScalarField<T>* field, uint64_t index)
  {
    if (_type == 0)
      return gradientMLS(grid, field, index);
  }

  template<typename T>
  Matrix<T> Approximator<T>::constructBMatrix(Grid<T>* grid,
    std::vector<uint64_t>* neighbors, uint64_t index, uint64_t order)
  {
    std::vector<std::vector<double> > B;
    for (uint64_t i = 0; i < neighbors->size(); i++)
    {
      uint64_t id = (*neighbors)[i];
      std::vector<double> temp = taylorMonomialExpansion(grid->getPoint(index),
            grid->getPoint(id), order);
      temp.push_back(1.0);
      B.push_back(temp);
    }
    return Matrix<T>("B",B);
  }

  template<typename T>
  std::vector<T> Approximator<T>::gradientMLS(Grid<T>* grid,
    ScalarField<T>* field, uint64_t index)
  {
    std::vector<T> result;
    //  TODO:: implement MLS here.
    //  First, find the nearest neighbors associated to the point specified by
    //  index.
    grid->queryNeighbors(_params.k);
    std::vector<uint64_t>* neighbors = grid->getNeighbors(index);
    //  Construct the matrix associated with the grid spacing
    Matrix<T> B = constructBMatrix(grid, neighbors, index, _params.n);
    return result;
  }
}
