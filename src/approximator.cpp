//  Approximator c++
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
