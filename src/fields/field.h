//------------------------------------------------------------------------------
//  field.h
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
#include <memory>

#include "grid.h"
#include "log.h"

//------------------------------------------------------------------------------
//  Forward declaration of Interpolator, Integrator and DiffEQ
//------------------------------------------------------------------------------
namespace ET
{
  //template<typename T> class ScalarField;
  template<typename T> class Interpolator;
  template<typename T> class DiffEQ;
  template<typename T> class Integrator;
}

//#include "scalarfield.h"
#include "interpolator.h"
#include "diffeq.h"
#include "integrator.h"

namespace ET
{
  //! \enum Field Types
  /*! An enum to collect the various types of grids.
   */
  enum class FieldType
  {
    DEFAULT,
    SCALAR,
    VECTOR,
    FRAME
  };
  //! Field Class
  /*! A Base class for various types of fields, such as ET::ScalarField,
   *  ET::VectorField and ET::FrameField.
   */
  template<typename T>
  class Field
  {
  public:
    //! Defualt Constructor
    /*! Default constructor for Field.
     */
    Field();
    //! Destructor
    /*! Destructor for Field.
     */
    ~Field();
    //! Constructor
    /*! constructor for Field that takes a name
     */
    Field(std::string t_name);
    //! Constructor
    /*! constructor for Field that takes a Logger
     */
    Field(std::shared_ptr<Log> t_log);
    //! Constructor
    /*! constructor for Field that takes a name and a Logger
     */
    Field(std::string t_name, std::shared_ptr<Log> t_log);
    //! Constructor
    /*! constructor for Field that takes a Grid
     */
    Field(std::shared_ptr<Grid<T>> t_grid);
    //! Constructor
    /*! constructor for Field that takes a name and a Grid
     */
    Field(std::string t_name, std::shared_ptr<Grid<T>> t_grid);
    //! Constructor
    /*! constructor for Field that takes a Interpolator
     */
    Field(std::shared_ptr<Interpolator<T>> t_interpolator);
    //! Constructor
    /*! constructor for Field that takes a name and a Interpolator
     */
    Field(std::string t_name, std::shared_ptr<Interpolator<T>> t_interpolator);
    //! Constructor
    /*! constructor for Field that takes a Grid and a logger
     */
    Field(std::shared_ptr<Grid<T>> t_grid, std::shared_ptr<Log> t_log);
    //! Constructor
    /*! constructor for Field that takes a name, Grid and a logger
     */
    Field(std::string, std::shared_ptr<Grid<T>> t_grid,
          std::shared_ptr<Log> t_log);

    //! Get Name
    /*! get name.  Get the name of the field.
     *  @return  m_name The name of the field.
     */
    std::string getName() const;
    //! Get Dim
    /*! get dimension.  Get the dimension of the field.
     *  @return  m_dim The dimension of the field.
     */
    size_t getDim() const;
    //! Get N
    /*! get number of elements.  Get the number of elements of the field.
     *  @return  m_N The number of elements of the field.
     */
    size_t getN() const;
    /*! Get Log.  Get the shared instance of the Log.
     *  @return A shared pointer to the Log.
     */
    std::shared_ptr<Log> getLog() const;
    /*! Get Grid.  Get the shared instance of the Grid.
     *  @return A shared pointer to the Grid.
     */
    std::shared_ptr<Grid<T>> getGrid() const;
    /*! Get Interpolator.  Get the shared instance of the Interpolator.
     *  @return A shared pointer to the Interpolator.
     */
    std::shared_ptr<Interpolator<T>> getInterpolator() const;
    /*! Get DiffEQ.  Get the shared instance of the DiffEQ.
     *  @return A shared pointer to the DiffEQ.
     */
    std::shared_ptr<DiffEQ<T>> getDiffEQ() const;
    /*! Get Integrator.  Get the shared instance of the Integrator.
     *  @return A shared pointer to the Integrator.
     */
    std::shared_ptr<Integrator<T>> getIntegrator() const;
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
    //! Get type
    /*! Get type.  Get the type of field.
     *  @return An FieldType enum
     */
    virtual enum FieldType getType() const;
    //! Set name
    /*! set name.  Sets the name of the Field.
        @param t_name a std::string for the name of the Field.
    */
    void setName(std::string t_name);
    //! Set Dim.
    /*! Set Dim.  Sets the dimension of the Field.
     *  @param t_dim The dimension to set the field to.
     */
    void setDim(size_t t_dim);
    //! Set N
    /*! Set N.  Sets the value of N for the field.
     *  @param t_N The number of points to set the Field to.
     */
    void setN(size_t t_N);
    //! Set Grid
    /*! set Grid.  Sets the shared pointer for the associated Grid.
     *  @param t_grid A shared pointer for a Grid instance.
     */
    void setGrid(std::shared_ptr<Grid<T>> t_grid);
    //! Set Log
    /*! set Log.  Sets the shared pointer for the associated Log.
     *  @param t_log A shared pointer for a Log instance.
     */
    void setLog(std::shared_ptr<Log> t_log);
    //! Set Interpolator
    /*! set Interpolator.  Sets the shared pointer for the
     *  associated Interpolator.
     *  @param t_interpolator A shared pointer for a Interpolator instance.
     */
    virtual void setInterpolator(std::shared_ptr<Interpolator<T>> t_interpolator);
    //! Set DiffEQ
    /*! set DiffEQ.  Sets the shared pointer for the associated DiffEQ.
     *  @param t_diffeq A shared pointer for a DiffEQ instance.
     */
    virtual void setDiffEQ(std::shared_ptr<DiffEQ<T>> t_diffeq);
    //! Set Integrator
    /*! set Integrator.  Sets the shared pointer for the associated Integrator.
     *  @param t_integrator A shared pointer for a Integrator instance.
     */
    virtual void setIntegrator(std::shared_ptr<Integrator<T>> t_integrator);
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
    //  Operators
    //! Assignment operator
    /*! Virtual function for the assignment operator which prevents slicing.
     */
    virtual Field<T>& operator=(const Field<T>& t_field) {
      assign(t_field);
      return *this;
    }
    //! Construct local field values
    /*! Function for generating vectors and matrices
     *  of local field values to use for interpolation.
     *  @param t_index The index of the point to construct around.
     *  @return Either vectors or matrices.
     */
    virtual Vector<T> constructLocalFieldValues(size_t t_index);
    //! Construct local field values
    /*! Function for generating vectors and matrices
     *  of local field values to use for interpolation.
     *  @param t_point The point to construct around.
     *  @param t_k The number of neighbors to use.
     *  @return Either vectors or matrices.
     */
    virtual Vector<T> constructLocalFieldValues(const std::vector<T>& t_point,
                                                size_t t_k);

  protected:
    /*! Name.  The name of the Interpolator.
     */
    std::string m_name {""};
    /*! Dim.  Dimension of the field.
     */
    size_t m_dim {0};
    /*! N.  Number of samples of the field.
     */
    size_t m_N {0};
    /*! Logger.  A shared instance of a Logger.
     */
    std::shared_ptr<Log> m_log;
    /*! Grid.  A shared instance of a Grid.
     */
    std::shared_ptr<Grid<T>> m_Grid;
    /*! Interpolator.  A shared instance of an Interpolator.
     */
    std::shared_ptr<Interpolator<T>> m_Interpolator;
    /*! DiffEQ.  A shared instance of a DiffEQ.
     */
    std::shared_ptr<DiffEQ<T>> m_DiffEQ;
    /*! Integrator.  A shared instance of an Integrator.
     */
    std::shared_ptr<Integrator<T>> m_Integrator;
    /*! Flag.  An in that denotes the type of message stored in info.
     */
    int m_flag;
    /*! Info.  A container for storing important information.
     */
    std::string m_info;
    /*! Field Type.  A const enum declaring the type of the Field.
     *  Defaulted to default.
     */
    const enum FieldType m_type {FieldType::DEFAULT};
    /*! Assign operation.  Copies the contents of t_field to this.
     */
    void assign(const Field<T>& t_field) {
      m_name = t_field.getName();
      m_dim = t_field.getDim();
      m_N = t_field.getN();
      m_log = t_field.getLog();
      m_Grid = t_field.getGrid();
      m_Interpolator = t_field.getInterpolator();
      m_DiffEQ = t_field.getDiffEQ();
      m_Integrator = t_field.getIntegrator();
      m_flag = t_field.getFlag();
      m_info = t_field.getInfo();
    }
  };

  template class Field<double>;

}
