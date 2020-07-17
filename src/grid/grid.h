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

#include "utils.h"
#include "matrix.h"
#include "kdtree.h"
#include "log.h"

#include <nanoflann.hpp>
#include "KDTreeVectorOfVectorsAdaptor.h"

namespace ET
{
  //! \class Grid Class
  /*! A Base class for various different types of grids.
   */
  template<typename T>
  class Grid
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
    //! Get Log
    /*! Get the shared logger for this Grid.
     *  @return m_log The shared instance of the logger.
     */
    std::shared_ptr<Log> getLog() const;
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
    //! Set Log
    /*! Set the shared logger for this Grid.
     *  @param t_log A std::shared_ptr<Log> instance of a logger.
     */
    void setLog(std::shared_ptr<Log> t_log);

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
    /*! Logger.  Shared instance of a logger.
     */
    std::shared_ptr<Log> m_log;

  };

  template class Grid<double>;

}
