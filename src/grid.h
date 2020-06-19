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
    Grid(uint64_t dim);
    Grid(std::string name, uint64_t dim);
    Grid(uint64_t dim, uint64_t N);
    Grid(std::string name, uint64_t dim, uint64_t N);

    //  Getters
    uint64_t getDim();
    uint64_t getN();
    std::vector<std::vector<T> >  getGrid();
    std::string getName();
    std::vector<std::vector<size_t> > getNeighbors();
    std::vector<std::vector<double> > getDistances();
    std::vector<std::vector<size_t> > getNeighborsRadius();
    std::vector<std::vector<double> > getDistancesRadius();
    std::vector<size_t>* getNeighbors(uint64_t index);

    //  Setters
    void setDim(uint64_t dim);
    void setN(uint64_t N);
    void setGrid(std::vector<std::vector<T> > grid);
    void setName(std::string name);

    //  Access operators for grid
    T& operator()(const uint64_t i, const uint64_t j);
    const T& operator()(const uint64_t i, const uint64_t j) const;
    //  Access operators for points
    std::vector<T>& operator()(const uint64_t i);
    const std::vector<T>& operator()(const uint64_t i) const;

    //  points and projections
    std::vector<T> getPoint(uint64_t i);
    std::vector<T> projection(uint64_t j);
    void setPoint(uint64_t i, std::vector<T> p);

    //  find nearest neighbors
    void queryNeighbors(uint64_t k);
    void queryRadius(double radius);

  private:
    uint64_t _dim;
    uint64_t _N;
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
