//------------------------------------------------------------------------------
//  Interpolator.h
//  The Entropic Trajectories Framework
//  -----------------------------------
//  Copyright (C) [2020] by [N. Carrara]
//  [ncarrara@albany.edu]

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
#include <iostream>
#include <map>
#include <memory>

#include "ugrid.h"
#include "geometry.h"
#include "params.h"
#include "utils.h"
#include "matrix.h"

//------------------------------------------------------------------------------
//  Forward declaration of ScalarField
//------------------------------------------------------------------------------
namespace ET
{
  template<typename T> class Field;
  //template<typename T> class ScalarField;
  //template<typename T> class VectorField;
}

#include "field.h"
//#include "scalarfield.h"
//#include "vectorfield.h"

namespace ET
{
  //! LSDriver enum
  /*! Enum for determining the type of least squares driver to use.
   *
   */
  enum class LSDriver
  {
    /*! Enum value: ET::LSDriver::xGELS.*/
    xGELS,
    /*! Enum value: ET::LSDriver::xGELSY.
     *  This uses a complete orthogonal transformation.
     */
    xGELSY,
    /*! Enum value: ET::LSDriver::xGELSD.
     *  This uses the SVD algorithm with divide and conquer.
     */
    xGELSD,
    /*! Enum value: ET::LSDriver::xGELSS.
     *  This uses standard SVD to solve the system.
     */
    xGELSS,
  };

  //! \class Interpolator Class
  /*! A Base class for interpolating on unstructured grids.
   *  Each instance has an associated shared pointer to a UGrid which
   *  gets passed to each derived class.
   */
  template<typename T>
  class Interpolator
  {
  public:
    //! Defualt Constructor
    /*! Default constructor for Interpolator.
     */
    Interpolator();
    //! Destructor
    /*! Destructor for Interpolator.
     */
    ~Interpolator();
    //! Constructor
    /*! constructor for Interpolator that takes a UGrid
     */
    Interpolator(std::shared_ptr<UGrid<T>> t_ugrid);
    //! Constructor
    /*! constructor for Interpolator that takes a Logger
     */
    Interpolator(std::shared_ptr<Log> t_log);
    //! Constructor
    /*! constructor for Interpolator that takes a UGrid and a logger
     */
    Interpolator(std::shared_ptr<UGrid<T>> t_ugrid,
                 std::shared_ptr<Log> t_log);

    /*! Get name.  Get the name of the Interpolator.
     *  @return The name of the Interpolator.
     */
    std::string getName() const;
    /*! Get UGrid.  Get the shared instance of the UGrid.
     *  @return A shared pointer to the UGrid.
     */
    std::shared_ptr<UGrid<T>> getUGrid() const;
    /*! Get Field.  Get the shared instance of the Field.
     *  @return A shared pointer to the Field.
     */
    std::shared_ptr<Field<T>> getField() const;
    /*! Get log.  Get the shared instance of the logger.
     *  @return A shared pointer to the logger
     */
    std::shared_ptr<Log> getLogger() const;
    //! Get flag
    /*! get flag.
     *  @return An int that classifies the type of information stored in
     *  info.
     */
    int getFlag() const;
    //! Get info
    /*! get info.
     *  @return The std::string info that contains relevant information.
     */
    std::string getInfo() const;
    //! Get LSDriver
    /*! get least squares driver.
     *  @return The LSDriver enum that is currently set to be used.
     */
    LSDriver getLSDriver() const;
    //! Set name
    /*! set name.  Sets the name of the Vector.
        @param t_name a std::string for the name of the Field.
    */
    void setName(std::string t_name);
    //! Set UGrid
    /*! set UGrid.  Sets the shared pointer for the associated UGrid.
     *  @param t_ugrid A shared pointer for a UGrid instance.
     */
    void setUGrid(std::shared_ptr<UGrid<T>> t_ugrid);
    //! Set Field
    /*! set Field.  Sets the shared pointer for the associated Field.
     *  @param t_field A shared pointer for a Field instance.
     */
    void setField(std::shared_ptr<Field<T>> t_field);
    //! Set Logger
    /*! set Logger.  Sets the shared pointer for the associated logger.
     *  @param t_log A shared pointer for a Log instance.
     */
    void setLogger(std::shared_ptr<Log> t_log);
    //! Set flag
    /*! set flag.  Sets the flag pertaining to info.
        @param t_flag an int that classifies the type of information stored
        in m_info.
    */
    void setFlag(int t_flag);
    //! Set info
    /*! set info.  Sets useful information pertaining to Field.
        @param t_info an std::string containing useful messages.
    */
    void setInfo(std::string t_info);
    //! Set LSDriver
    /*! set least squares driver.  Sets the type of the least squares
     *  method to use.
     *  @param t_type A string denoting the type of LS driver to use.
     */
    void setLSDriver(std::string type);

    //  Derivative functions that must be overloaded in derived classes.

    //! Derivative
    /*! derivative.  Derivative for a point in the UGrid given by index,
     *  of degree t_degree.
     *  @return The nth-derivative at the point given
     *  by the index.
     */
    virtual Vector<T> derivative(const size_t t_index,
                                 const size_t t_degree);
    //! Derivative
    /*! derivative.  Derivative for a point in the UGrid given by index,
     *  of degree t_degree and in the direction t_direction.
     *  @return The nth-derivative in the lth-direction at the point given
     *  by the index.
     */
    virtual T derivative(const size_t t_index,
                         const size_t t_degree,
                         const size_t t_direction);
   //! Derivative
   /*! derivative.  Derivative for an arbitrary point,
    *  of degree t_degree.
    *  @return The nth-derivative at the point given
    *  by the index.
    */
   virtual Vector<T> derivative(const std::vector<T>& point,
                                const size_t t_degree);
   //! Derivative
   /*! derivative.  Derivative for an arbitrary point,
    *  of degree t_degree and in the direction t_direction.
    *  @return The nth-derivative in the lth-direction at the point given
    *  by the index.
    */
   virtual T derivative(const std::vector<T>& point,
                        const size_t t_degree,
                        const size_t t_direction);
    //! xGELSx
    /*! xGELSx.  Helper function for running the least squares
     *  algorithm for the system Ax = y.
     *  @param t_A The Matrix A.
     *  @param t_x The Vector x.
     *  @return The solution Vector y.
     */
    Vector<T> xGELSx(Matrix<T> A, Vector<T> x);
    //! xGELSx
    /*! xGELSx.  Helper function for running the least squares
     *  algorithm for the system AX = Y.
     *  @param t_A The Matrix A.
     *  @param t_X The Matrix X
     *  @return The solution Matrix Y.
     */
    Matrix<T> xGELSx(Matrix<T> A, Matrix<T> X);
    //! xMLSx
    /*! xMLSx.  Helper function for running the moving least squares
     *  algorithm for the system Ax = y.
     *  @param t_A The Matrix A.
     *  @param t_x The Vector x.
     *  @return The solution Vector y.
     */
    Vector<T> xMLSx(Matrix<T> A, Vector<T> x);
    //! xMLSx
    /*! xMLSx.  Helper function for running the moving least squares
     *  algorithm for the system AX = Y.
     *  @param t_A The Matrix A.
     *  @param t_X The Matrix X
     *  @return The solution Matrix Y.
     */
    Matrix<T> xMLSx(Matrix<T> A, Matrix<T> X);

    //--------------------------------------------------------------------------
    //  various functions
    //--------------------------------------------------------------------------
    const std::string summary();
    //--------------------------------------------------------------------------

  protected:
    /*! Name.  The name of the Interpolator. */
    std::string m_name {""};
    /*! Logger.  A shared instance of a Logger.*/
    std::shared_ptr<Log> m_log;
    /*! UGrid.  A shared instance of a UGrid.*/
    std::shared_ptr<UGrid<T>> m_ugrid;
    /*! Field.  A shared instance of a Field.*/
    std::shared_ptr<Field<T>> m_field;
    /*! LSDriver.  An enum that denotes the type of least squares
     *  method to use.
     */
    enum LSDriver m_lsdriver {LSDriver::xGELS};
    /*! Flag.  An in that denotes the type of message stored in info.*/
    int m_flag {0};
    /*! Info.  A container for storing important information.*/
    std::string m_info {""};
  };

  template class Interpolator<double>;
}
