//  Class for an unstructured grid
#pragma once

#include <vector>
#include <string>
#include "matrix.h"

namespace ET
{
  template<typename T>
  class Grid
  {
  public:
    Grid();
    ~Grid();
    Grid(unsigned int dim);
    Grid(std::string name, unsigned int dim);
    Grid(unsigned int dim, unsigned int N);
    Grid(std::string name, unsigned int dim, unsigned int N);

    //  Getters
    unsigned int get_dim();
    unsigned int get_N();
    Matrix<T> get_grid();
    std::string get_name();

    //  Setters
    void set_dim(unsigned int dim);
    void set_N(unsigned int N);
    void set_grid(Matrix<T> grid);
    void set_name(std::string name);

    //  Access operators for grid
    T& operator()(const unsigned int i, const unsigned int j);
    const T& operator()(const unsigned int i, const unsigned int j) const;

    //  points and projections
    std::vector<T> get_point(unsigned int i);
    std::vector<T> projection(unsigned int i);
    void set_point(unsigned int i, std::vector<T> point);

  private:
    unsigned int _dim;
    unsigned int _N;
    Matrix<T> _grid;
    std::string _name;
    std::vector<std::string> _coords;
  };

}
