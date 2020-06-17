//  Class for an unstructured grid
#pragma once

#include <vector>
#include <string>
#include "utils.h"
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
    unsigned int getDim();
    unsigned int getN();
    std::vector<std::vector<T> >  getGrid();
    std::string getName();
    std::vector<std::vector<size_t> > getNeighbors();
    std::vector<std::vector<double> > getDistances();
    std::vector<std::vector<size_t> > getNeighborsRadius();
    std::vector<std::vector<double> > getDistancesRadius();

    //  Setters
    void setDim(unsigned int dim);
    void setN(unsigned int N);
    void setGrid(std::vector<std::vector<T> > grid);
    void setName(std::string name);

    //  Access operators for grid
    T& operator()(const unsigned int i, const unsigned int j);
    const T& operator()(const unsigned int i, const unsigned int j) const;
    //  Access operators for points
    std::vector<T>& operator()(const unsigned int i);
    const std::vector<T>& operator()(const unsigned int i) const;

    //  points and projections
    std::vector<T> getPoint(unsigned int i);
    std::vector<T> projection(unsigned int j);
    void setPoint(unsigned int i, std::vector<T> p);

    //  find nearest neighbors
    void queryNeighbors(unsigned int k);
    void queryRadius(double radius);

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
    std::vector<std::vector<size_t> > _neighbors_radius;
    std::vector<std::vector<double> > _distances_radius;
    //KDTreeVectorOfVectorsAdaptor<std::vector<std::vector<T> >, T> _kdt;
  };

  template class Grid<double>;

}
