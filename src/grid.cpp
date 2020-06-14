//  Grid cpp
#include "grid.h"

namespace ET
{
  template<typename T>
  Grid<T>::Grid()
  {

  }
  template<typename T>
  Grid<T>::~Grid()
  {

  }
  template<typename T>
  Grid<T>::Grid(unsigned int dim) : _dim(dim)
  {

  }
  template<typename T>
  Grid<T>::Grid(std::string name, unsigned int dim) : _dim(dim), _name(name)
  {

  }
  template<typename T>
  Grid<T>::Grid(unsigned int dim, unsigned int N) : _dim(dim), _N(N)
  {

  }
  template<typename T>
  Grid<T>::Grid(std::string name, unsigned int dim, unsigned int N)
  : _dim(dim), _N(N), _name(name)
  {

  }
  template<typename T>
  unsigned int Grid<T>::get_dim()
  {
    return _dim;
  }
  template<typename T>
  unsigned int Grid<T>::get_N()
  {
    return _N;
  }
  template<typename T>
  std::vector<Point<T> > Grid<T>::get_grid()
  {
    return _grid;
  }
  template<typename T>
  std::string Grid<T>::get_name()
  {
    return _name;
  }

  //  Setters
  template<typename T>
  void Grid<T>::set_dim(unsigned int dim)
  {
    _dim = dim;
  }
  template<typename T>
  void Grid<T>::set_N(unsigned int N)
  {
    _N = N;
  }
  template<typename T>
  void Grid<T>::set_grid(Matrix<T> grid)
  {
    _grid = grid;
  }
  template<typename T>
  void Grid<T>::set_name(std::string name)
  {
    _name = name;
  }

  //  Access operators for grid
  template<typename T>
  T& Grid<T>::operator()(const unsigned int i, const unsigned int j)
  {
    return _grid(i,j);
  }
  template<typename T>
  const T& Grid<T>::operator()(const unsigned int i, const unsigned int j) const
  {
    return _grid(i,j);
  }

  //  points and projections
  template<typename T>
  std::vector<Point<T> > Grid<T>::get_point(unsigned int i)
  {
    return _grid.get_row(i);
  }
  template<typename T>
  std::vector<T> Grid<T>::projection(unsigned int i)
  {
    return _grid.get_col(i);
  }
  template<typename T>
  void Grid<T>::set_point(unsigned int i, std::vector<Point<T> > point)
  {
    _grid.set_row(i,point);
  }

}
