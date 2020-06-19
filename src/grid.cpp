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
  Grid<T>::Grid(uint64_t dim) : _dim(dim)
  {

  }
  template<typename T>
  Grid<T>::Grid(std::string name, uint64_t dim) : _dim(dim), _name(name)
  {

  }
  template<typename T>
  Grid<T>::Grid(uint64_t dim, uint64_t N) : _dim(dim), _N(N)
  {

  }
  template<typename T>
  Grid<T>::Grid(std::string name, uint64_t dim, uint64_t N)
  : _dim(dim), _N(N), _name(name)
  {

  }
  template<typename T>
  uint64_t Grid<T>::getDim()
  {
    return _dim;
  }
  template<typename T>
  uint64_t Grid<T>::getN()
  {
    return _N;
  }
  template<typename T>
  std::vector<std::vector<T> > Grid<T>::getGrid()
  {
    return _grid;
  }
  template<typename T>
  std::string Grid<T>::getName()
  {
    return _name;
  }
	template<typename T>
	std::vector<std::vector<size_t> > Grid<T>::getNeighbors()
	{
		return _neighbors;
	}
	template<typename T>
	std::vector<std::vector<double> > Grid<T>::getDistances()
	{
		return _distances;
	}
	template<typename T>
	std::vector<std::vector<size_t> > Grid<T>::getNeighborsRadius()
	{
		return _neighbors_radius;
	}
	template<typename T>
	std::vector<std::vector<double> > Grid<T>::getDistancesRadius()
	{
		return _distances_radius;
	}
	template<typename T>
	std::vector<size_t>* Grid<T>::getNeighbors(uint64_t index)
	{
		return &_neighbors[index];
	}
  //  Setters
  template<typename T>
  void Grid<T>::setDim(uint64_t dim)
  {
    _dim = dim;
  }
  template<typename T>
  void Grid<T>::setN(uint64_t N)
  {
    _N = N;
  }
  template<typename T>
  void Grid<T>::setGrid(std::vector<std::vector<T> > grid)
  {
    _grid = grid;
		_N = grid.size();
  }
  template<typename T>
  void Grid<T>::setName(std::string name)
  {
    _name = name;
  }

  //  Access operators for grid
  template<typename T>
  T& Grid<T>::operator()(const uint64_t i, const uint64_t j)
  {
    return _grid[i][j];
  }
  template<typename T>
  const T& Grid<T>::operator()(const uint64_t i, const uint64_t j) const
  {
    return _grid[i][j];
  }
	template<typename T>
  std::vector<T>& Grid<T>::operator()(const uint64_t i)
  {
    return _grid[i];
  }
  template<typename T>
  const std::vector<T>& Grid<T>::operator()(const uint64_t i) const
  {
    return _grid[i];
  }

  //  points and projections
  template<typename T>
  std::vector<T> Grid<T>::getPoint(uint64_t i)
  {
    return _grid[i];
  }
  template<typename T>
  std::vector<T> Grid<T>::projection(uint64_t j)
  {
		std::vector<T> result(_N);
		for (uint64_t i = 0; i < _N; i++)
		{
			result[i] = _grid[i][j];
		}
    return result;
  }
  template<typename T>
  void Grid<T>::setPoint(uint64_t i, std::vector<T> point)
  {
    _grid[i] = point;
  }

  template<typename T>
  void Grid<T>::queryNeighbors(uint64_t k)
  {
    _neighbors.resize(_N);
		_distances.resize(_N);
    //  generate kdtree
    KDTreeVectorOfVectorsAdaptor<std::vector<std::vector<T> >, T> kdt(_dim, _grid, 16);
    kdt.index->buildIndex();

    const size_t num_results = k;
    for (uint64_t i = 0; i < _N; i++)
    {
      std::vector<size_t> ret_indexes(num_results);
      std::vector<double> out_dists_sqr(num_results);

      nanoflann::KNNResultSet<double> resultSet(num_results);
      resultSet.init(&ret_indexes[0], &out_dists_sqr[0] );
			kdt.index->findNeighbors(resultSet, &_grid[i][0], nanoflann::SearchParams(10));
			_neighbors[i] = std::move(ret_indexes);
			_distances[i] = std::move(out_dists_sqr);
    }
  }

	template<typename T>
  void Grid<T>::queryRadius(double radius)
  {
    _neighbors_radius.resize(_N);
		_distances_radius.resize(_N);
		//	the algorithm looks for points that satisfy the squared
		//	distance rather than the square root.
		radius *= radius;
    //  generate kdtree
    KDTreeVectorOfVectorsAdaptor<std::vector<std::vector<T> >, T> kdt(_dim, _grid, 16);
    kdt.index->buildIndex();

    for (uint64_t i = 0; i < _N; i++)
    {
      std::vector<std::pair<size_t,double> > ret_matches;

			kdt.index->radiusSearch(&_grid[i][0], radius, ret_matches, nanoflann::SearchParams(10));
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
