//------------------------------------------------------------------------------
//  grid.h
//  The Entropic Trajectories Framework
//  -----------------------------------
//  Copyright (C) [2020] by [N. Carrara]
//  [ncarrara@albany.edu]
//
//  Permission to use, copy, modify, and/or distribute this software for any
//  purpose with or without fee t_is hereby granted.
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
#include <limits>

#include "utilities.h"
#include "matrix.h"
#include "kdtree.h"
#include "log.h"

#include <nanoflann.hpp>
#include "KDTreeVectorOfVectorsAdaptor.h"

namespace ET
{

  //! \enum Grid Types
  /*! An enum to collect the various types of grids.
   */
  enum class GridType
  {
    DEFAULT,
    STRUCTURED,
    UNSTRUCTURED,
    MESHLESS,
  };
  //! \enum Boundary Types
  /*! An enum for a boundary type for a point.
   *  The boundary type specifies implicit boundaries
   *  for each degree of freedom.
   */
  enum class BoundaryType
  {
    UNBOUNDED,
    OPEN,
    CLOSED,
  };
  //! \class BoundaryPair
  /*! A struct for holding boundary information for
   *  each degree of freedom.
   */
  template<typename T>
  struct BoundaryPair
  {
    //! Left boundarytype
    enum BoundaryType m_left_type = {BoundaryType::UNBOUNDED};
    //! Right boundarytype
    enum BoundaryType m_right_type = {BoundaryType::UNBOUNDED};
    //! Left bounary value
    T m_left = -std::numeric_limits<T>::infinity();
    //! Right boundary value
    T m_right = std::numeric_limits<T>::infinity();
    //! Constructor
    BoundaryPair() {}
    //! Constructor with values (assuming both sides are closed)
    BoundaryPair(const T& t_left, const T& t_right)
    : m_left(t_left), m_right(t_right),
      m_left_type(BoundaryType::CLOSED), m_right_type(BoundaryType::CLOSED)
    {
    }
    //! Constructor with values and types
    BoundaryPair(const T& t_left, const T& t_right,
                 const enum BoundaryType t_left_type,
                 const enum BoundaryType t_right_type)
    : m_left(t_left), m_right(t_right),
      m_left_type(t_left_type), m_right_type(t_right_type)
    {
    }
  };
  //! \enum
  //! \class Grid Class
  /*! A Base class for various different types of grids.
   */
  template<typename T>
  class Grid : public std::enable_shared_from_this<Grid<T>>
  {
  public:
    //  Constructors
    //! Defualt Constructor
    /*! Default constructor for Grid.
     */
    Grid();
    //! Destructor
    /*! Destructor for Grid.
     */
    ~Grid();
    //! Constructor
    /*! constructor for Grid that takes a Logger
     *  @param t_log A shared logger instance.
     */
    Grid(std::shared_ptr<Log> t_log);
    //! Constructor
    /*! constructor for Grid that takes a name
     *  @param t_name An std::string specifying this objects name.
     */
    Grid(std::string t_name);
    //! Constructor
    /*! constructor for Grid that takes a name and a Logger
     *  @param t_name An std::string specifying this objects name.
     *  @param t_log A shared logger instance.
     */
    Grid(std::string t_name, std::shared_ptr<Log> t_log);
    //! Constructor
    /*! constructor for Grid that takes a dimension.
     *  @param t_dim A size_t object which specifies the dimension.
     */
    Grid(size_t t_dim);
    //! Constructor
    /*! constructor for Grid that takes a dimension and a logger.
     *  @param t_dim A size_t object which specifies the dimension.
     *  @param t_log A shared logger instance.
     */
    Grid(size_t t_dim, std::shared_ptr<Log> t_log);
    //! Constructor
    /*! constructor for Grid that takes a name and a dimension.
     *  @param t_name An std::string specifying this objects name.
     *  @param t_dim A size_t object which specifies the dimension.
     */
    Grid(std::string t_name, size_t t_dim);
    //! Constructor
    /*! constructor for Grid that takes a name, a dimension and a logger.
     *  @param t_name An std::string specifying this objects name.
     *  @param t_dim A size_t object which specifies the dimension.
     *  @param t_log A shared logger instance.
     */
    Grid(std::string t_name, size_t t_dim, std::shared_ptr<Log> t_log);
    //! Constructor
    /*! constructor for Grid that takes a dimension and a number of elements.
     *  @param t_dim A size_t object which specifies the dimension.
     *  @param t_N A size_t object which specifies the number of elements.
     */
    Grid(size_t t_dim, size_t t_N);
    //! Constructor
    /*! constructor for Grid that takes a dimension, a number of elements
     *  and a logger.
     *  @param t_dim A size_t object which specifies the dimension.
     *  @param t_log A shared logger instance.
     */
    Grid(size_t t_dim, size_t t_N, std::shared_ptr<Log> t_log);
    //! Constructor
    /*! constructor for Grid that takes a name, a dimension and
     *  a number of elements.
     *  @param t_name An std::string specifying this objects name.
     *  @param t_dim A size_t object which specifies the dimension.
     */
    Grid(std::string t_name, size_t t_dim, size_t t_N);
    //! Constructor
    /*! constructor for Grid that takes a name, a dimension,
     *  a number of elementsand a logger.
     *  @param t_name An std::string specifying this objects name.
     *  @param t_dim A size_t object which specifies the dimension.
     *  @param t_log A shared logger instance.
     */
    Grid(std::string t_name, size_t t_dim, size_t t_N,
         std::shared_ptr<Log> t_log);
    //! Constructor
    /*! constructor for Grid that takes a grid array
    *  @param t_grid A std::vector<std::vector<T>> array of points.
    *  @param move_grid A bool that determines whether the t_grid object
    *  should be copied or moved to m_grid.
    */
    Grid(std::vector<std::vector<T>> t_grid, bool move_grid=false);
    //! Constructor
    /*! constructor for Grid that takes a grid array and a logger
    *  @param t_grid A std::vector<std::vector<T>> array of points.
    *  @param t_log A shared logger instance.
    *  @param move_grid A bool that determines whether the t_grid object
    *  should be copied or moved to m_grid.
    */
    Grid(std::vector<std::vector<T>> t_grid, std::shared_ptr<Log> t_log,
         bool move_grid=false);
    //! Constructor
    /*! constructor for Grid that takes a name and a grid array
    *  @param t_name An std::string specifying this objects name.
    *  @param t_grid A std::vector<std::vector<T>> array of points.
    *  @param move_grid A bool that determines whether the t_grid object
    *  should be copied or moved to m_grid.
    */
    Grid(std::string t_name, std::vector<std::vector<T>> t_grid,
         bool move_grid=false);
    //! Constructor
    /*! constructor for Grid that takes a name, grid array and a logger
    *  @param t_name An std::string specifying this objects name.
    *  @param t_grid A std::vector<std::vector<T>> array of points.
    *  @param t_log A shared logger instance.
    *  @param move_grid A bool that determines whether the t_grid object
    *  should be copied or moved to m_grid.
    */
    Grid(std::string t_name, std::vector<std::vector<T>> t_grid,
         std::shared_ptr<Log> t_log, bool move_grid=false);
    //   Getters and Setters
    //! Get name.
    /*! Get the name of the Grid.
     *  @return m_name The name of the Grid.
     */
    std::string getName() const;
    //! Get dim
    /*! Get the dimension of the Grid.
     *  @return m_dim The dimension of the grid.
     */
    size_t getDim() const;
    //! Get N
    /*! Get the number of elements in the Grid.
     *  @return m_N The number of elements in the grid.
     */
    size_t getN() const;
    //! Get Grid
    /*! Get the array containing the points in the grid.
     *  @return m_grid The array containing the points in the grid.
     */
    std::vector<std::vector<T>> getGrid() const;
    //! Get Coords
    /*! Get the list of coordinate labels.
     *  @return m_coords The list of coordinate labels.
     */
    std::vector<std::string> getCoords() const;
    //! Get Log
    /*! Get the shared logger for this Grid.
     *  @return m_log The shared instance of the logger.
     */
    std::shared_ptr<Log> getLog() const;
    //! Get type
    /*! Get the type of the grid.  This will be useful for other classes
     *  which only have a shared pointer to a generic Grid<T>, but do not
     *  otherwise know which functions are available.
     */
    enum GridType getType() const;
    //! Get BoundaryPairs
    /*! Get the std::vector of boundary pairs.
     *  @return The std::vector of boundary pairs
     */
    std::vector<BoundaryPair<T>> getBoundaryPairs() const;
    //! Get BoundaryPair
    /*! Get a specific boundary pair.
     *  @param t_i The index of the boundary pair.
     *  @return A BoundaryPair struct.
     */
    BoundaryPair<T> getBoundaryPair(size_t t_i) const;
    //! Get point
    /*! Get a particular point in the grid.  This method is identical to
     *  the operator()(const size_t t_i)
     *  @param t_i The index of the desired point
     *  @return A std::vector<T> of the point.
     */
    std::vector<T> getPoint(const size_t t_i) const;
    //! Set name.
    /*! Set the name of the Grid.
     *  @param t_name A std::string specifying the name of the Grid.
     */
    void setName(const std::string t_name);
    //! Set dim
    /*! Set the dimension of the Grid.
     *  @param t_dim A size_t object specifying the dimension of the Grid.
     */
    void setDim(const size_t t_dim);
    //! Set N
    /*! Set the number of elements in the Grid.
     *  @param t_N A size_t object specifying the number of
     *  elements in the Grid.
     */
    void setN(const size_t t_N);
    //! Set Grid
    /*! Set the array containing the points in the grid.
     *  @param t_grid A std::vector<std::vector<T>> of the points.
     */
    void setGrid(const std::vector<std::vector<T>> t_grid);
    //! Set Coords
    /*! Set the coordinate labels for the grid.
     *  @param t_coords A list of coordinate labels for the grid.
     */
    void setCoords(const std::vector<std::string> t_coords);
    //! Set Log
    /*! Set the shared logger for this Grid.
     *  @param t_log A std::shared_ptr<Log> instance of a logger.
     */
    void setLog(std::shared_ptr<Log> t_log);
    //! Set BoundaryPairs
    /*! Set the std::vector of boundary pairs.
     *  @param The std::vector of boundary pairs
     */
    void setBoundaryPairs(std::vector<BoundaryPair<T>> t_boundary);
    //! Set BoundaryPair
    /*! Set a specific boundary pair.
     *  @param t_i The index of the boundary pair.
     *  @param A BoundaryPair struct.
     */
    void setBoundaryPair(size_t t_i, BoundaryPair<T> t_boundary);
    //! Move data
    /*! Reassign ownership of the elements of an std::vector<std::vector<T>>
     *  to the m_grid data member.
     *  @param t_grid A std::vector<std::vector<T>> of data.
     */
    void moveGrid(const std::vector<std::vector<T>> t_grid);

    //  Operator overloads
    /*!
        @param t_Grid A const Grid<T>& reference.
        @return A copy of the Grid
    *///--------------------------------------------------------------------------
    Grid<T>& operator=(const Grid<T>& t_Grid);
    /*! Is equal.  Determines whether two grids are equivalent.
        @param t_Grid A const Grid<T>& reference.
        @return A boolean quantifying whether this Grid and Grid
        are equivalent.
    */
    bool operator==(const Grid<T>& t_Grid) const;
    /*! Is not equal.  Determines whether two grids are not equal.
        @param t_Grid A const Grid<T>& reference.
        @return A boolean quantifying whether this Grid and Grid
        are not equivalent.
    */
    bool operator!=(const Grid<T>& t_Grid) const;
    /*! Minus.  Constructs a copy of the Grid \f$\vec{v}\f$ multiplied by
        minus one, \f$-\vec{v}\f$.
        @return A copy of this Grid multiplied by \f$-1\f$.
    */
    Grid<T> operator-() const;
    /*! Sum.  Adds two grids together.
        @param t_Grid A const Grid<T>& reference.
        @return The sum of this Grid and Grid.
    */
    Grid<T> operator+(const Grid<T>& t_Grid) const;
    /*! Sum equals.  Adds a Grid to the current one.
        @param t_Grid A const Grid<T>& reference.
        @return This Grid plus Grid
    */
    Grid<T>& operator+=(const Grid<T>& t_Grid);
    /*! Difference.  Subtracts two grids.
        @param t_Grid A const Grid<T>& reference.
        @return The difference of this Grid and Grid.
    */
    Grid<T> operator-(const Grid<T>& t_Grid) const;
    /*! Difference equals.  Subtracts a Grid from the current one.
        @param t_Grid A const Grid<T>& reference.
        @return This Grid minus Grid.
    */
    Grid<T>& operator-=(const Grid<T>& t_Grid);
    //  Access operators
    /*! Access operator.  An access operator for changing entries of the Grid.
        @param t_i A const size_t& reference for the index of
        the desired element.
        @param t_j A const size_t& reference for the index of
        the desired dimension.
        @return The (i,j)-element of the array.
    */
    T& operator()(const size_t& t_i, const size_t& t_j);
    /*! Access operator.  An access operator for retrieving coefficients
        of the Grid.
        @param t_i A const size_t& reference for the index of
        the desired element.
        @param t_j A const size_t& reference for the index of
        the desired dimension.
        @return The (i,j)-element of the array.
    */
    const T& operator()(const size_t& t_i, const size_t& t_j) const;
    //  Access operators
    /*! Access operator.  An access operator for changing entries of the Grid.
        @param t_i A const size_t& reference for the index of
        the desired element.
        @return The ith-point in the Grid.
    */
    std::vector<T>& operator()(const size_t& t_i);
    /*! Access operator.  An access operator for retrieving coefficients
        of the Grid.
        @param t_i A const size_t& reference for the index of
        the desired element.
        @return The ith-point in the Grid.
    */
    const std::vector<T>& operator()(const size_t& t_i) const;

    //  Special functions
    //! Projection
    /*! Proj. Get a projection along a particular axis.
     *  @param t_axis The axis to project along.
     *  @return An std::vector<T> of the coordinates along that axis.
     */
    std::vector<T> proj(const size_t t_axis);
    //! Projection
    /*! Proj.  Get a projection along a set of axes.
     *  @param t_axes The axes to project along.
     *  @return An std::vector<std::vector<T>> of the coordinates along
     *  the projection.
     */
    std::vector<std::vector<T>> proj(const std::vector<size_t> t_axes);

    //  Virtual functions for subclasses
    //! get KDTree
    virtual std::shared_ptr<KDTree<T>> getKDTree() const { return std::make_shared<KDTree<T>>(); }
    //! set KDtree
    virtual void setKDTree(std::shared_ptr<KDTree<T>> t_kdtree) { return; }
    //! Get shared ptr.
    /*! Returns a shared_ptr of this instance.
     */
    std::shared_ptr<Grid<T>> getptr() {
      return std::shared_ptr<Grid<T>>(this);
    }

  protected:
    /*! Name.  Name of the Grid.  Defaulted to empty string.
     */
    std::string m_name {""};
    /*! Dimension.  Dimension of the grid.  Defaulted to 0.
     */
    size_t m_dim {0};
    /*! N.  Number of elements in the grid.  Default to 0.
     */
    size_t m_N {0};
    /*! grid.  Array containing the points of the grid.
     *  Defaulted to an single point whose value is zero.
     */
    std::vector<std::vector<T>> m_grid {{{0}}};
    /*! Coordinates.  List of values labeling the coordinates.
     *  Defaulted to list containing a single empty string.
     */
    std::vector<std::string> m_coords {{""}};
    /*! Logger.  Shared instance of a logger.
     */
    std::shared_ptr<Log> m_log {std::make_shared<Log>()};
    /*! Grid Type.  A const enum declaring the type of the grid.
     *  Defaulted to default.
     */
    const enum GridType m_type {GridType::DEFAULT};
    /*! Vector of boundary pairs for each dof.
     */
    std::vector<BoundaryPair<T>> m_boundary;
    /*! Logging system name generator.
     */
    virtual std::string NAME() const {
      return "Grid:[" + std::to_string(m_id.id) + "]:" + m_name + ":";
    }
    /*! Unique ID for each instance.
     */
    UniqueID m_id;
  public:
    //  Inheriting iterators from std::vector
    using iterator = typename std::vector<std::vector<T>>::iterator;
    using const_iterator = typename std::vector<std::vector<T>>::const_iterator;

    iterator begin()              { return m_grid.begin(); }
    iterator end()                { return m_grid.end(); }
    const_iterator cbegin() const { return m_grid.cbegin(); }
    const_iterator cend() const   { return m_grid.cend(); }
  };

  template struct BoundaryPair<double>;
  template class Grid<double>;

}
