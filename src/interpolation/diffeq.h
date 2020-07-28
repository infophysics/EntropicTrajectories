//------------------------------------------------------------------------------
//  diffeq.h
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

#include<vector>
#include<string>
#include <iostream>
#include <map>
#include <memory>

#include "grid.h"
#include "ugrid.h"
#include "log.h"
#include "utilities.h"

// ------------------------------------------------------------------------------
//  Forward declaration of Interpolator, Integrator and ScalarField
// ------------------------------------------------------------------------------
namespace ET
{
  template<typename T> class Field;
  template<typename T> class Interpolator;
  //template<typename T> class ScalarField;
  template<typename T> class Integrator;
}
#include "field.h"
#include "interpolator.h"
//#include "scalarfield.h"
#include "integrator.h"

namespace ET
{
  //! \enum DiffEQ Type Enum
  /*! An anum for quantifying the type of differential eq.
   */
  enum class DiffEQType
  {
    DEFAULT,
    FIRST_ORDER_ODE,
    SECOND_ORDER_ODE,
    PDE,
  };

  //! \class DiffEQ
  /*! A base class for differential equations
   */
  template<typename T>
  class DiffEQ : public std::enable_shared_from_this<DiffEQ<T>>
  {
  public:
    //! Default Constructor
    /*
     */
    DiffEQ();
    //! Destructor
    /*!
     */
    ~DiffEQ();
    //! Constructor with a shared log instance
    /*!
     */
    DiffEQ(std::shared_ptr<Log> t_log);
    //! Constructor with a Grid and a log
    /*!
     */
    DiffEQ(std::shared_ptr<Grid<T>> t_Grid,
           std::shared_ptr<Log> t_log);

    //  Getters and Setters
    //! Get Name
    /*! Get the name of the DiffEQ.
     */
    std::string getName() const;
    //! Get dimension.
    /*! Get the dimension of the DiffEQ.
     */
    size_t getDim() const;
    //! Get Grid
    /*! Get the shared_ptr for the Grid.
     */
    std::shared_ptr<Grid<T>> getGrid() const;
    //! Get Log
    /*! Get the shared_ptr for the Log.
     */
    std::shared_ptr<Log> getLog() const;
    //! Get Type
    /*! Get the type of the DiffEQ.
     */
    enum DiffEQType getType() const;
    //! Set Name
    /*! Set the name of the DiffEQ.
     */
    void setName(const std::string t_name);
    //! Set Dimension.
    /*! Set the dimension of the DiffEQ.
     */
    void setDim(const size_t t_dim);
    //! Set Grid
    /*! Set the Grid of the DiffEQ.
     */
    void setGrid(const std::shared_ptr<Grid<T>> t_Grid);
    //! Set Field
    /*! Set the Field of the DiffEQ.
     */
    void setField(const std::shared_ptr<Field<T>> t_Field);
    //! Set the Log
    /*! Set the Log of the DiffEQ.
     */
    void setLog(const std::shared_ptr<Log> t_log);

    //  Overloaded functions
    //! Evaluation at a point
    /*! Evaluate the diffEQ at a point and time.
     */
    virtual T evaluate(const std::vector<T>& t_point, double t_time);
    //! Evaluation at many points
    /*! Evaluate the diffEQ at many points and a time
     */
    virtual std::vector<T> evaluate(const std::vector<std::vector<T>>& t_points,
                                    double t_time);
    //! Get shared ptr.
    /*! Returns a shared_ptr of this instance.
     */
    std::shared_ptr<DiffEQ<T>> getptr() {
      return std::shared_ptr<DiffEQ<T>>(this);
    }


  protected:
    /*! Name.  Name of the DiffEQ.  Defaulted to empty string.
     */
    std::string m_name {""};
    /*! Dimension.  Dimension of the DiffEQ.  Defaulted to 0.
     */
    size_t m_dim {0};
    /*! Grid.  A shared_ptr for the underlying Grid.
     */
    std::shared_ptr<Grid<T>> m_Grid;
    /*! Logger.  Shared instance of a logger.
     */
    std::shared_ptr<Log> m_log {std::make_shared<Log>()};
    /* Type of the DiffEQ
     */
    const enum DiffEQType m_type {DiffEQType::DEFAULT};
    /*! Logging system name generator.
     */
    virtual std::string NAME() const {
      return "DiffEQ:[" + std::to_string(m_id.id) + "]:" + m_name + ":";
    }
    /*! Unique ID for each instance.
     */
    UniqueID m_id;

  };

  //! \class FirstOrderODE
  /*! A base class for first order ODE's
   */
  template<typename T>
  class FirstOrderODE : public DiffEQ<T>
  {
  public:
    //! Default Constructor
    /*!
     */
    FirstOrderODE();
    //! Destructor
    /*!
     */
    ~FirstOrderODE();
    //! Constructor with a shared log instance
    /*!
     */
    FirstOrderODE(std::shared_ptr<Log> log);
    //! Constructor with a Grid and a log
    /*!
     */
    FirstOrderODE(std::shared_ptr<Grid<T>> uGrid, std::shared_ptr<Log> log);
    //  Overloaded functions
    //! Evaluation at a point
    /*! Evaluate the diffEQ at a point and time.
     */
    virtual T evaluate(const std::vector<T>& t_point, double t_time);
    //! Evaluation at many points
    /*! Evaluate the diffEQ at many points and a time
     */
    virtual std::vector<T> evaluate(const std::vector<std::vector<T>>& t_points,
                                    double t_time);

  protected:
    /* Type of the DiffEQ
     */
    const enum DiffEQType m_type {DiffEQType::FIRST_ORDER_ODE};
    /*! Logging system name generator.
     */
    virtual std::string NAME() const {
      return "FirstOrderODE:[" + std::to_string(this->m_id.id) + "]:" + this->m_name + ":";
    }
  };

  //! \class SecondOrderODE
  /*! A base class for second order ODE's
   */
  template<typename T>
  class SecondOrderODE : public DiffEQ<T>
  {
  public:
    //! Default Constructor
    /*!
     */
    SecondOrderODE();
    //! Destructor
    /*!
     */
    ~SecondOrderODE();
    //! Constructor with a shared log instance
    /*!
     */
    SecondOrderODE(std::shared_ptr<Log> log);
    //! Constructor with a Grid and a log
    /*!
     */
    SecondOrderODE(std::shared_ptr<Grid<T>> uGrid, std::shared_ptr<Log> log);
    //  Overloaded functions
    //! Evaluation at a point
    /*! Evaluate the diffEQ at a point and time.
     */
    virtual T evaluate(const std::vector<T>& t_point, double t_time);
    //! Evaluation at many points
    /*! Evaluate the diffEQ at many points and a time
     */
    virtual std::vector<T> evaluate(const std::vector<std::vector<T>>& t_points,
                                    double t_time);

  protected:
    /* Type of the DiffEQ
     */
    const enum DiffEQType m_type {DiffEQType::SECOND_ORDER_ODE};
    /*! Logging system name generator.
     */
    virtual std::string NAME() const {
      return "SecondOrderODE:[" + std::to_string(this->m_id.id) + "]:" + this->m_name + ":";
    }
  };

  //! \class PDE
  /*! A base class for partial differential EQ's
   */
  template<typename T>
  class PDE : public DiffEQ<T>
  {
  public:
    //! Default Constructor
    /*!
     */
    PDE();
    //! Destructor
    /*!
     */
    ~PDE();
    //! Constructor with a shared log instance
    /*!
     */
    PDE(std::shared_ptr<Log> log);
    //! Constructor with a Grid and a log
    /*!
     */
    PDE(std::shared_ptr<Grid<T>> uGrid, std::shared_ptr<Log> log);
    //  Overloaded functions
    //! Evaluation at a point
    /*! Evaluate the diffEQ at a point and time.
     */
    virtual T evaluate(const std::vector<T>& t_point, double t_time);
    //! Evaluation at many points
    /*! Evaluate the diffEQ at many points and a time
     */
    virtual std::vector<T> evaluate(const std::vector<std::vector<T>>& t_points,
                                    double t_time);

  protected:
    /* Type of the DiffEQ
     */
    const enum DiffEQType m_type {DiffEQType::PDE};
    /*! Logging system name generator.
     */
    virtual std::string NAME() const {
      return "PDE:[" + std::to_string(this->m_id.id) + "]:" + this->m_name + ":";
    }
  };


  template class DiffEQ<double>;
  template class FirstOrderODE<double>;
  template class SecondOrderODE<double>;
  template class PDE<double>;

}
