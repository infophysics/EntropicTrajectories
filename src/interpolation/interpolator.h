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

#include "grid.h"
#include "geometry.h"
#include "params.h"
#include "utilities.h"
#include "matrix.h"
#include "log.h"

//------------------------------------------------------------------------------
//  Forward declaration of ScalarField
//------------------------------------------------------------------------------
namespace ET
{
  template<typename T> class Field;
  //template<typename T> class ScalarField;
  template<typename T> class DiffEQ;
  template<typename T> class Integrator;
  //template<typename T> class VectorField;
}

#include "field.h"
//#include "scalarfield.h"
#include "diffeq.h"
#include "integrator.h"
//#include "vectorField.h"

namespace ET
{
  //! \enum Interpolator Type
  /*! Enum for determining the type of interpolation scheme.
   */
  enum class InterpolatorType
  {
    /*! Enum value: ET::InterpolatorType::DEFAULT.*/
    DEFAULT,
    /*! Enum value: ET::InterpolatorType::LTE.*/
    LTE,
    /*! Enum value: ET::InterpolatorType::RBF.*/
    RBF,
  };
  //! \enum LSDriver enum
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
  //! \enum Solver Type enum
  /*! Determines whether one uses LS, MLS or WMLS methods for
   *  various algorithms.
   */
  enum class SolverType
  {
    /*! Enum value: ET::SolverType::LS.*/
    LS,
    /*! Enum value: ET::SolverType::MLS.*/
    MLS,
    /*! Enum value: ET::SolverType::WMLS.*/
    WMLS,
  };
  //! \enum Weight matrix type
  /*! Enum for determining the type of weight matrix to use
   *  in the WMLS method.
   */
  enum class WeightMatrixType
  {
    /*! Enum value: ET::WeightMatrixType::GAUSSIAN.*/
    GAUSSIAN,
  };

  extern std::map<std::string, InterpolatorType> InterpolatorTypeMap;
  extern std::map<InterpolatorType, std::string> InterpolatorTypeNameMap;
  extern std::map<std::string, LSDriver> LSDriverMap;
  extern std::map<LSDriver, std::string> LSDriverNameMap;
  extern std::map<std::string, SolverType> SolverTypeMap;
  extern std::map<SolverType, std::string> SolverTypeNameMap;
  extern std::map<std::string, WeightMatrixType> WeightMatrixTypeMap;
  extern std::map<WeightMatrixType, std::string> WeightMatrixTypeNameMap;

  //! \class Interpolator Class
  /*! A Base class for interpolating on unstructured Grids.
   *  Each instance has an associated shared pointer to a Grid which
   *  gets passed to each derived class.
   */
  template<typename T>
  class Interpolator : public std::enable_shared_from_this<Interpolator<T>>
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
    /*! constructor for Interpolator that takes a name
     */
    Interpolator(std::string t_name);
    //! Constructor
    /*! constructor for Interpolator that takes a Grid
     */
    Interpolator(std::shared_ptr<Grid<T>> t_Grid);
    //! Constructor
    /*! constructor for Interpolator that takes a name and a Grid
     */
    Interpolator(std::string t_name, std::shared_ptr<Grid<T>> t_Grid);
    //! Constructor
    /*! constructor for Interpolator that takes a Log
     */
    Interpolator(std::shared_ptr<Log> t_log);
    //! Constructor
    /*! constructor for Interpolator that takes a name and a Log
     */
    Interpolator(std::string t_name, std::shared_ptr<Log> t_log);
    //! Constructor
    /*! constructor for Interpolator that takes a Grid and a logger
     */
    Interpolator(std::shared_ptr<Grid<T>> t_Grid,
                 std::shared_ptr<Log> t_log);
    //! Constructor
    /*! constructor for Interpolator that takes a name, Grid and a logger
    */
    Interpolator(std::string t_name, std::shared_ptr<Grid<T>> t_Grid,
                 std::shared_ptr<Log> t_log);

    /*! Get name.  Get the name of the Interpolator.
     *  @return The name of the Interpolator.
     */
    std::string getName() const;
    /*! Get Grid.  Get the shared instance of the Grid.
     *  @return A shared pointer to the Grid.
     */
    std::shared_ptr<Grid<T>> getGrid() const;
    /*! Get Field.  Get the shared instance of the Field.
     *  @return A shared pointer to the Field.
     */
    std::shared_ptr<Field<T>> getField() const;
    /*! Get log.  Get the shared instance of the logger.
     *  @return A shared pointer to the logger
     */
    std::shared_ptr<Log> getLog() const;
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
    //! Get SolverType
    /*! get the solver type.
     *  @return The SolverType enum that is currently set.
     */
    SolverType getSolverType() const;
    //! Get LSDriver
    /*! get least squares driver.
     *  @return The LSDriver enum that is currently set to be used.
     */
    LSDriver getLSDriver() const;
    //! Get InterpolatorType
    /*! get least squares driver.
     *  @return The InterpolatorType enum that is currently set to be used.
     */
    InterpolatorType getInterpolatorType() const;
    //! Set name
    /*! set name.  Sets the name of the Vector.
        @param t_name a std::string for the name of the Field.
    */
    void setName(std::string t_name);
    //! Set Grid
    /*! set Grid.  Sets the shared pointer for the associated Grid.
     *  @param t_Grid A shared pointer for a Grid instance.
     */
    void setGrid(std::shared_ptr<Grid<T>> t_Grid);
    //! Set Field
    /*! set Field.  Sets the shared pointer for the associated Field.
     *  @param t_Field A shared pointer for a Field instance.
     */
    void setField(std::shared_ptr<Field<T>> t_Field);
    //! Set Log
    /*! set Log.  Sets the shared pointer for the associated logger.
     *  @param t_log A shared pointer for a Log instance.
     */
    void setLog(std::shared_ptr<Log> t_log);
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
    void setLSDriver(std::string t_type);
    //! Set SolverType
    /*! set the solver type.  Sets the solver type to use,
     *  either LS, MLS or WMLS.
     *  @param t_type A string denoting the solver type.
     */
    void setSolverType(std::string t_type);

    //  Derivative functions that must be overloaded in derived classes.
    //! Derivative
    /*! derivative.  Derivative for a point in the Grid given by index,
     *  of degree t_degree.
     *  @return The nth-derivative at the point given
     *  by the index.
     */
    virtual Vector<T> derivative(const size_t t_index,
                                 const size_t t_degree);
    //! Derivative
    /*! derivative.  Derivative for a point in the Grid given by index,
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
    //! Field Derivative
    /*! field derivative.  Derivative for the entire field
     *  of degree t_degree.
     *  @return The nth-derivative of the entire field
     */
    virtual std::vector<Vector<T>> fieldDerivative(const size_t t_degree);
    //! Field Derivative
    /*! field derivative.  Derivative for the entire field
     *  of degree t_degree and in the direction t_direction.
     *  @return The nth-derivative in the lth-direction of the entire field
     */
    virtual std::vector<T> fieldDerivative(const size_t t_degree,
                                           const size_t t_direction);
    //! xLinearSolvex
    /*! xLinearSolvex.  Helper function for running either,
     *  LS, MLS or WMLS methods.
     *  @param t_A The Matrix A.
     *  @param t_x The Vector x.
     *  @param t_index Index of the point in the grid (only for WMLS)
     *  @return The solution Vector y.
     */
    Vector<T> xLinearSolvex(Matrix<T> A, Vector<T> x,
                            size_t t_index = (size_t)-1);
    //! xLinearSolvex
    /*! xLinearSolvex.  Helper function for running either,
     *  LS, MLS or WMLS methods.
     *  @param t_A The Matrix A.
     *  @param t_X The Matrix X.
     *  @param t_index Index of the point in the grid (only for WMLS)
     *  @return The solution Matrix Y.
     */
    Matrix<T> xLinearSolvex(Matrix<T> A, Matrix<T> x,
                            size_t t_index = (size_t)-1);
    //! xLinearSolvex
    /*! xLinearSolvex.  Helper function for running either,
     *  LS, MLS or WMLS methods.
     *  @param t_A The Matrix A.
     *  @param t_x The Vector x.
     *  @param t_point The point to expand around (only for WMLS)
     *  @return The solution Vector y.
     */
    Vector<T> xLinearSolvex(Matrix<T> A, Vector<T> x,
                            std::vector<T> t_point = {});
    //! xLinearSolvex
    /*! xLinearSolvex.  Helper function for running either,
     *  LS, MLS or WMLS methods.
     *  @param t_A The Matrix A.
     *  @param t_X The Matrix X.
     *  @param t_point The point to expand around (only for WMLS)
     *  @return The solution Matrix Y.
     */
    Matrix<T> xLinearSolvex(Matrix<T> A, Matrix<T> x,
                            std::vector<T> t_point = {});
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
    //! xWMLSx
    /*! xWMLSx.  Helper function for running the weighted moving least squares
     *  algorithm for the system Ax = y.
     *  @param t_point The point to expand around.
     *  @param t_A The Matrix A.
     *  @param t_x The Vector x.
     *  @param t_index The index to expand around.
     *  @return The solution Vector y.
     */
    Vector<T> xWMLSx(Matrix<T> A, Vector<T> x, size_t t_index);
    //! xWMLSx
    /*! xWMLSx.  Helper function for running the weighted moving least squares
     *  algorithm for the system AX = Y.
     *  @param t_A The Matrix A.
     *  @param t_X The Matrix X.
     *  @param t_index The index to expand around.
     *  @return The solution Matrix Y.
     */
    Matrix<T> xWMLSx(Matrix<T> A, Matrix<T> X, size_t t_index);
    //! xWMLSx
    /*! xWMLSx.  Helper function for running the weighted moving least squares
     *  algorithm for the system Ax = y.
     *  @param t_A The Matrix A.
     *  @param t_x The Vector x.
     *  @param t_point The point to expand around.
     *  @return The solution Vector y.
     */
    Vector<T> xWMLSx(Matrix<T> A, Vector<T> x, std::vector<T> t_point);
    //! xWMLSx
    /*! xWMLSx.  Helper function for running the weighted moving least squares
     *  algorithm for the system AX = Y.
     *  @param t_A The Matrix A.
     *  @param t_X The Matrix X.
     *  @param t_point The point to expand around.
     *  @return The solution Matrix Y.
     */
    Matrix<T> xWMLSx(Matrix<T> A, Matrix<T> X, std::vector<T> t_point);

  protected:
    /*! Name.  The name of the Interpolator.
     */
    std::string m_name {""};
    /*! Log.  A shared instance of a Log.
     */
    std::shared_ptr<Log> m_log;
    /*! Grid.  A shared instance of a Grid.
     */
    std::shared_ptr<Grid<T>> m_Grid;
    /*! Field.  A shared instance of a Field.
     */
    std::shared_ptr<Field<T>> m_Field;
    /*! LSDriver.  An enum that denotes the type of least squares
     *  method to use.
     */
    enum LSDriver m_lsdriver {LSDriver::xGELS};
    /*! Solver Type.  An enum that denotes the type of solver to use.
     */
    enum SolverType m_solvertype {SolverType::LS};
    /*! Flag.  An in that denotes the type of message stored in info.
     */
    int m_flag {0};
    /*! Info.  A container for storing important information.
     */
    std::string m_info {""};
    /*! Interpolator Type.  A const identifier of the type of interpolator.
     */
    const enum InterpolatorType m_interpolatortype {InterpolatorType::DEFAULT};
    /*! Logging system name generator.
     */
    virtual std::string NAME() const {
      return "Interpolator:[" + std::to_string(m_id.id) + "]:" + m_name + ":";
    }
    /*! Unique ID for each instance.
     */
    UniqueID m_id;
  };

  template class Interpolator<double>;
}
