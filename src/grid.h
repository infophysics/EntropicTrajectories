//  Class for an unstructured grid
#pragma once

#include <vector>
#include <string>
#include "matrix.h"
#include <nanoflann.hpp>
#include "KDTreeVectorOfVectorsAdaptor.h"


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
    std::vector<std::vector<T> >  get_grid();
    std::string get_name();
    std::vector<std::vector<size_t> > get_neighbors();
    std::vector<std::vector<double> > get_distances();

    //  Setters
    void set_dim(unsigned int dim);
    void set_N(unsigned int N);
    void set_grid(std::vector<std::vector<T> > grid);
    void set_name(std::string name);

    //  Access operators for grid
    T& operator()(const unsigned int i, const unsigned int j);
    const T& operator()(const unsigned int i, const unsigned int j) const;

    //  points and projections
    std::vector<std::vector<T> > get_point(unsigned int i);
    std::vector<T> projection(unsigned int i);
    void set_point(unsigned int i, std::vector<std::vector<T> > p);

    //  find nearest neighbors
    void find_neighbors(unsigned int k);

  private:
    unsigned int _dim;
    unsigned int _N;
    //  Unstructured grid
    std::vector<std::vector<T> > _grid;
    std::string _name;
    std::vector<std::string> _coords;
    //  associated kdTree neighbors
    std::vector<std::vector<size_t> > _neighbors;
    std::vector<std::vector<double> > _distances;

    //KDTreeVectorOfVectorsAdaptor<std::vector<std::vector<T> >, T> _kdt;
  };

}
