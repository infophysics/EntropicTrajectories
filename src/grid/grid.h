//------------------------------------------------------------------------------
//  grid.h
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
