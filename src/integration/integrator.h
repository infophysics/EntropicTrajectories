//------------------------------------------------------------------------------
//  integrator.h
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
#include "ugrid.h"
#include "params.h"
#include "utilities.h"
#include "matrix.h"
//------------------------------------------------------------------------------
//  Forward declaration of ScalarField
//------------------------------------------------------------------------------
namespace ET
{
  template<typename T> class Field;
  //template<typename T> class ScalarField;
  template<typename T> class Interpolator;
  template<typename T> class DiffEQ;
  //template<typename T> class VectorField;
}
#include "field.h"
//#include "scalarfield.h"
#include "diffeq.h"
#include "interpolator.h"
//#include "vectorfield.h"

namespace ET
{
  //! \enum Integrator type enum
  /*!
   */
  enum class IntegratorType
  {
    DEFAULT,
    ODE,
    PDE,
  };

  //! \enum ODE Integration Scheme Enum
  /*!
   */
  enum class ODEIntegrationScheme
  {
    EULER,
    IMPROVED_EULER,
    RUNGE_KUTTA_4,
  };

  //! \enum PDE Integration Scheme Enum
  /*!
   */
  enum class PDEIntegrationScheme
  {
    FDM,
  };

  //! \class Integrator Base class.
  /*!
   */
  template<typename T>
  class Integrator
  {
  public:
    //! Default Constructor
    /*
     */
    Integrator();
    //! Destructor
    /*!
     */
    ~Integrator();
    //! Constructor with a shared log instance
    /*!
     */
    Integrator(std::shared_ptr<Log> t_log);
    //! Constructor with a Grid and a log
    /*!
     */
    Integrator(std::shared_ptr<Grid<T>> t_Grid, std::shared_ptr<Log> t_log);
    //  Getters and Setters
    //! Get Name
    /*! Get the name of the Integrator.
     */
    std::string getName() const;
    //! Get dimension.
    /*! Get the dimension of the Integrator.
     */
    size_t getDim() const;
    //! Get Grid
    /*! Get the shared_ptr for the Grid.
     */
    std::shared_ptr<Grid<T>> getGrid() const;
    //! Get DiffEQ
    /*! Get the shared_ptr for the DiffEQ.
     */
    std::shared_ptr<DiffEQ<T>> getDiffEQ() const;
    //! Get Log
    /*! Get the shared_ptr for the Log.
     */
    std::shared_ptr<Log> getLog() const;
    //! Get Type
    /*! Get the type of the Integrator.
     */
    enum IntegratorType getType() const;
    //! Set Name
    /*! Set the name of the Integrator.
     */
    void setName(const std::string t_name);
    //! Set Dimension.
    /*! Set the dimension of the Integrator.
     */
    void setDim(const size_t t_dim);
    //! Set Grid
    /*! Set the Grid of the Integrator.
     */
    void setGrid(const std::shared_ptr<Grid<T>> t_Grid);
    //! Set DiffEQ
    /*! Set the DiffEQ of the Integrator.
     */
    void setDiffEQ(const std::shared_ptr<DiffEQ<T>> t_DiffEQ);
    //! Set the Log
    /*! Set the Log of the Integrator.
     */
    void setLog(const std::shared_ptr<Log> t_log);
    //! Get shared ptr.
    /*! Returns a shared_ptr of this instance.
     */
    std::shared_ptr<Integrator<T>> getptr() {
      return std::shared_ptr<Integrator<T>>(this);
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
    /*! DiffEQ.  Shared pointer to the DiffEQ.
     */
    std::shared_ptr<DiffEQ<T>> m_DiffEQ;
    /*! Type.  The type of the integrator
     */
    const enum IntegratorType m_type {IntegratorType::DEFAULT};
    /*! Logging system name generator.
     */
    virtual std::string NAME() const {
      return "Integrator:[" + std::to_string(m_id.id) + "]:" + m_name + ":";
    }
    /*! Unique ID for each instance.
     */
    UniqueID m_id;
  };

  //! \class ODEIntegrator Base class.
  /*!
   */
  template<typename T>
  class ODEIntegrator : Integrator<T>
  {
  public:
    //! Default Constructor
    /*
     */
    ODEIntegrator();
    //! Destructor
    /*!
     */
    ~ODEIntegrator();
    //! Constructor with a shared log instance
    /*!
     */
    ODEIntegrator(std::shared_ptr<Log> t_log);
    //! Constructor with a Grid and a log
    /*!
     */
    ODEIntegrator(std::shared_ptr<Grid<T>> t_Grid, std::shared_ptr<Log> t_log);
    //! Get ODE Scheme
    /*! Get the scheme of the Integrator.
     */
    enum ODEIntegrationScheme getScheme() const;
    //! Set ODE Scheme
    /*! Set the scheme for the ODE Integrator
     */
    void setScheme(const enum ODEIntegrationScheme t_scheme);

  protected:
    /*! Type.  The type of the integrator
     */
    const enum IntegratorType m_type {IntegratorType::ODE};
    /*! ODE Scheme.  The ODE scheme to use
     */
    enum ODEIntegrationScheme m_scheme {ODEIntegrationScheme::EULER};
    /*! Logging system name generator.
     */
    virtual std::string NAME() const {
      return "ODEIntegrator:[" + std::to_string(this->m_id.id)
             + "]:" + this->m_name + ":";
    }
  };

  //! \class PDEIntegrator Base class.
  /*!
   */
  template<typename T>
  class PDEIntegrator : Integrator<T>
  {
  public:
    //! Default Constructor
    /*
     */
    PDEIntegrator();
    //! Destructor
    /*!
     */
    ~PDEIntegrator();
    //! Constructor with a shared log instance
    /*!
     */
    PDEIntegrator(std::shared_ptr<Log> t_log);
    //! Constructor with a Grid and a log
    /*!
     */
    PDEIntegrator(std::shared_ptr<Grid<T>> t_Grid, std::shared_ptr<Log> t_log);
    //! Get PDE Scheme
    /*! Get the scheme of the Integrator.
     */
    enum PDEIntegrationScheme getScheme() const;
    //! Set PDE Scheme
    /*! Set the scheme for the PDE Integrator
     */
    void setScheme(const enum PDEIntegrationScheme t_scheme);

  protected:
    /*! Type.  The type of the integrator
     */
    const enum IntegratorType m_type {IntegratorType::PDE};
    /*! PDE Scheme.  The PDE scheme to use
     */
    enum PDEIntegrationScheme m_scheme {PDEIntegrationScheme::FDM};
    /*! Logging system name generator.
     */
    virtual std::string NAME() const {
      return "PDEIntegrator:[" + std::to_string(this->m_id.id)
             + "]:" + this->m_name + ":";
    }
  };

  template class Integrator<double>;
}
