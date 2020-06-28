//------------------------------------------------------------------------------
//  ugrid.cpp
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
#include "ugrid.h"

namespace ET
{
	//----------------------------------------------------------------------------
  //  UGrid constructors
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Default constructor
  //    sets name = " ", and _dim, _N = 0
  //----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::UGrid() : _dim(0), _N(0)
  {
  }
  template<typename T>
  UGrid<T>::~UGrid()
  {
  }
  template<typename T>
  UGrid<T>::UGrid(uint64_t dim) : _dim(dim)
  {
  }
  template<typename T>
  UGrid<T>::UGrid(std::string name, uint64_t dim) : _dim(dim), _name(name)
  {
  }
  template<typename T>
  UGrid<T>::UGrid(uint64_t dim, uint64_t N) : _dim(dim), _N(N)
  {
  }
  template<typename T>
  UGrid<T>::UGrid(std::string name, uint64_t dim, uint64_t N)
  : _dim(dim), _N(N), _name(name)
  {
  }
	//	One dimensional constructor
	template<typename T>
	UGrid<T>::UGrid(std::vector<T> ugrid)
	{
		_N = ugrid.size();
		_ugrid.resize(_N);
		_dim = 1;
		for (uint32_t i = 0; i < _N; i++)
		{
			std::vector<T> temp = {ugrid[i]};
			_ugrid[i] = temp;
		}
	}
	//	n-dimensional constructor
	template<typename T>
	UGrid<T>::UGrid(std::vector<std::vector<T> > ugrid)
	{
		_N = ugrid.size();
		_dim = ugrid[0].size();
		_ugrid = ugrid;
	}
  template<typename T>
  uint64_t UGrid<T>::getDim()
  {
    return _dim;
  }
  template<typename T>
  uint64_t UGrid<T>::getN()
  {
    return _N;
  }
  template<typename T>
  std::vector<std::vector<T> > UGrid<T>::getUGrid()
  {
    return _ugrid;
  }
  template<typename T>
  std::string UGrid<T>::getName()
  {
    return _name;
  }
	template<typename T>
	std::vector<std::vector<size_t> > UGrid<T>::getNeighbors()
	{
		return _neighbors;
	}
	template<typename T>
	std::vector<std::vector<double> > UGrid<T>::getDistances()
	{
		return _distances;
	}
	template<typename T>
	std::vector<std::vector<size_t> > UGrid<T>::getNeighborsRadius()
	{
		return _neighbors_radius;
	}
	template<typename T>
	std::vector<std::vector<double> > UGrid<T>::getDistancesRadius()
	{
		return _distances_radius;
	}
	template<typename T>
	std::vector<size_t>* UGrid<T>::getNeighbors(uint64_t index)
	{
		return &_neighbors[index];
	}
  //  Setters
  template<typename T>
  void UGrid<T>::setDim(uint64_t dim)
  {
    _dim = dim;
  }
  template<typename T>
  void UGrid<T>::setN(uint64_t N)
  {
    _N = N;
  }
  template<typename T>
  void UGrid<T>::setUGrid(std::vector<std::vector<T> > ugrid)
  {
    _ugrid = ugrid;
		_N = ugrid.size();
  }
  template<typename T>
  void UGrid<T>::setName(std::string name)
  {
    _name = name;
  }

  //  Access operators for ugrid
  template<typename T>
  T& UGrid<T>::operator()(const uint64_t i, const uint64_t j)
  {
    return _ugrid[i][j];
  }
  template<typename T>
  const T& UGrid<T>::operator()(const uint64_t i, const uint64_t j) const
  {
    return _ugrid[i][j];
  }
	template<typename T>
  std::vector<T>& UGrid<T>::operator()(const uint64_t i)
  {
    return _ugrid[i];
  }
  template<typename T>
  const std::vector<T>& UGrid<T>::operator()(const uint64_t i) const
  {
    return _ugrid[i];
  }

  //  points and projections
  template<typename T>
  std::vector<T> UGrid<T>::getPoint(uint64_t i)
  {
    return _ugrid[i];
  }
  template<typename T>
  std::vector<T> UGrid<T>::projection(uint64_t j)
  {
		std::vector<T> result(_N);
		for (uint64_t i = 0; i < _N; i++)
		{
			result[i] = _ugrid[i][j];
		}
    return result;
  }
  template<typename T>
  void UGrid<T>::setPoint(uint64_t i, std::vector<T> point)
  {
    _ugrid[i] = point;
  }

  template<typename T>
  void UGrid<T>::queryNeighbors(uint64_t k)
  {
    _neighbors.resize(_N);
		_distances.resize(_N);
    //  generate kdtree
    KDTreeVectorOfVectorsAdaptor<std::vector<std::vector<T> >, T> kdt(_dim, _ugrid, 16);
    kdt.index->buildIndex();

    const size_t num_results = k;
    for (uint64_t i = 0; i < _N; i++)
    {
      std::vector<size_t> ret_indexes(num_results);
      std::vector<double> out_dists_sqr(num_results);

      nanoflann::KNNResultSet<double> resultSet(num_results);
      resultSet.init(&ret_indexes[0], &out_dists_sqr[0] );
			kdt.index->findNeighbors(resultSet, &_ugrid[i][0], nanoflann::SearchParams(10));
			_neighbors[i] = std::move(ret_indexes);
			_distances[i] = std::move(out_dists_sqr);
    }
  }

	template<typename T>
  void UGrid<T>::queryRadius(double radius)
  {
    _neighbors_radius.resize(_N);
		_distances_radius.resize(_N);
		//	the algorithm looks for points that satisfy the squared
		//	distance rather than the square root.
		radius *= radius;
    //  generate kdtree
    KDTreeVectorOfVectorsAdaptor<std::vector<std::vector<T> >, T> kdt(_dim, _ugrid, 16);
    kdt.index->buildIndex();

    for (uint64_t i = 0; i < _N; i++)
    {
      std::vector<std::pair<size_t,double> > ret_matches;

			kdt.index->radiusSearch(&_ugrid[i][0], radius, ret_matches, nanoflann::SearchParams(10));
			std::vector<size_t> indices(ret_matches.size());
			std::vector<double> distances(ret_matches.size());
			for (uint64_t j = 0; j < ret_matches.size(); j++)
			{
				indices[j] = ret_matches[j].first;
				distances[j] = ret_matches[j].second;
			}
			_neighbors_radius[i] = std::move(indices);
			_distances_radius[i] = std::move(distances);
    }
  }

}
