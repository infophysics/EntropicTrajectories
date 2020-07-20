//------------------------------------------------------------------------------
//  ugrid.h
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

#include "grid.h"
#include "utils.h"
#include "matrix.h"
#include "kdtree.h"
#include "log.h"

#include <nanoflann.hpp>
#include "KDTreeVectorOfVectorsAdaptor.h"

namespace ET
{
  //----------------------------------------------------------------------------
  //  Class for unstructured grids
  //----------------------------------------------------------------------------
  template<typename T>
  class UGrid : public Grid<T>
  {
  public:
    //  Constructors
    //! Defualt Constructor
    /*! Default constructor for Grid.
     */
    UGrid();
    //! Destructor
    /*! Destructor for UGrid.
     */
    ~UGrid();
    //! Constructor
    /*! constructor for UGrid that takes a Logger
     *  @param t_log A shared logger instance.
     */
    UGrid(std::shared_ptr<Log> t_log);
    //! Constructor
    /*! constructor for UGrid that takes a name
     *  @param t_name An std::string specifying this objects name.
     */
    UGrid(std::string t_name);
    //! Constructor
    /*! constructor for UGrid that takes a name and a Logger
     *  @param t_name An std::string specifying this objects name.
     *  @param t_log A shared logger instance.
     */
    UGrid(std::string t_name, std::shared_ptr<Log> t_log);
    //! Constructor
    /*! constructor for UGrid that takes a dimension.
     *  @param t_dim A size_t object which specifies the dimension.
     */
    UGrid(size_t t_dim);
    //! Constructor
    /*! constructor for UGrid that takes a dimension and a logger.
     *  @param t_dim A size_t object which specifies the dimension.
     *  @param t_log A shared logger instance.
     */
    UGrid(size_t t_dim, std::shared_ptr<Log> t_log);
    //! Constructor
    /*! constructor for UGrid that takes a name and a dimension.
     *  @param t_name An std::string specifying this objects name.
     *  @param t_dim A size_t object which specifies the dimension.
     */
    UGrid(std::string t_name, size_t t_dim);
    //! Constructor
    /*! constructor for UGrid that takes a name, a dimension and a logger.
     *  @param t_name An std::string specifying this objects name.
     *  @param t_dim A size_t object which specifies the dimension.
     *  @param t_log A shared logger instance.
     */
    UGrid(std::string t_name, size_t t_dim, std::shared_ptr<Log> t_log);
    //! Constructor
    /*! constructor for UGrid that takes a dimension and a number of elements.
     *  @param t_dim A size_t object which specifies the dimension.
     *  @param t_N A size_t object which specifies the number of elements.
     */
    UGrid(size_t t_dim, size_t t_N);
    //! Constructor
    /*! constructor for UGrid that takes a dimension, a number of elements
     *  and a logger.
     *  @param t_dim A size_t object which specifies the dimension.
     *  @param t_log A shared logger instance.
     */
    UGrid(size_t t_dim, size_t t_N, std::shared_ptr<Log> t_log);
    //! Constructor
    /*! constructor for UGrid that takes a name, a dimension and
     *  a number of elements.
     *  @param t_name An std::string specifying this objects name.
     *  @param t_dim A size_t object which specifies the dimension.
     */
    UGrid(std::string t_name, size_t t_dim, size_t t_N);
    //! Constructor
    /*! constructor for UGrid that takes a name, a dimension,
     *  a number of elementsand a logger.
     *  @param t_name An std::string specifying this objects name.
     *  @param t_dim A size_t object which specifies the dimension.
     *  @param t_log A shared logger instance.
     */
    UGrid(std::string t_name, size_t t_dim, size_t t_N,
          std::shared_ptr<Log> t_log);
    //! Constructor
    /*! constructor for Grid that takes a grid array
    *  @param t_grid A std::vector<std::vector<T>> array of points.
    */
    UGrid(std::vector<std::vector<T>> t_grid);
    //! Constructor
    /*! constructor for Grid that takes a grid array and a logger
    *  @param t_grid A std::vector<std::vector<T>> array of points.
    *  @param t_log A shared logger instance.
    */
    UGrid(std::vector<std::vector<T>> t_grid, std::shared_ptr<Log> t_log);
    //! Constructor
    /*! constructor for Grid that takes a name and a grid array
    *  @param t_name An std::string specifying this objects name.
    *  @param t_grid A std::vector<std::vector<T>> array of points.
    */
    UGrid(std::string t_name, std::vector<std::vector<T>> t_grid);
    //! Constructor
    /*! constructor for Grid that takes a name, grid array and a logger
    *  @param t_name An std::string specifying this objects name.
    *  @param t_grid A std::vector<std::vector<T>> array of points.
    *  @param t_log A shared logger instance.
    */
    UGrid(std::string t_name, std::vector<std::vector<T>> t_grid,
          std::shared_ptr<Log> t_log);

    //  Getters and Setters
    //! Get KDTree.
    /*! Get the KDTree.
     *  @return A KDTree<T> object.
     */
    KDTree<T> getKDTree() const;
    //! Set KDTree.
    /*! Set the KDTree.
     *  @param t_kdtree A KDTree<T> object.
     */
    void setKDTree(KDTree<T> t_kdtree);

    //  Various functions
    bool checkConsistency();

  private:
    /*! KDTree for the grid.*/
    KDTree<T> m_kdtree;
    /*! Grid Type.  The grid type for this object is UNSTRUCTURED.*/
    const enum GridType t_type {GridType::UNSTRUCTURED};
  };

  template class UGrid<double>;

}
