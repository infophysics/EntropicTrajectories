//------------------------------------------------------------------------------
//  vectorfield.cpp
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
#include "vectorfield.h"


namespace ET
{
  //----------------------------------------------------------------------------
  //  VectorField constructors
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Default constructor
  //    sets name = "default", and _dim, _N = 0
  //----------------------------------------------------------------------------
  template<typename T>
  VectorField<T>::VectorField() : _dim(0), _N(0), _name("default")
  {
    //##########################################################################
    _log = std::make_shared<Log>();
    _log->init("ET:VectorField:default", ".logs/vectorfield_default.txt");
    _log->TRACE("Vector Field 'default' created at location "
                + getMem(*this));
    //##########################################################################
  }
  //----------------------------------------------------------------------------
  //  Default destructor
  //----------------------------------------------------------------------------
  template<typename T>
  VectorField<T>::~VectorField()
  {
    //##########################################################################
    _log->TRACE("Vector Field '" + _name
                + "' destroyed at location " + getMem(*this));
    //##########################################################################
  }
  //----------------------------------------------------------------------------
  //  Various constructors taking in arguments for
  //		_dim, _name, _N, _ugrid and _log.
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //	Constructor with UGrid
  //----------------------------------------------------------------------------
  template<typename T>
  VectorField<T>::VectorField(std::shared_ptr<UGrid<T>> ugrid)
  : _ugrid(ugrid), _name("default")
  {
    _N = _ugrid->getN();
    _dim = _ugrid->getDim();
    _approx = std::make_shared<Interpolator<T>>();
    //##########################################################################
    _log = std::make_shared<Log>();
    _log->init("ET:VectorField:default", ".logs/vectorfield_default.txt");
    _log->TRACE("Vector Field 'default' created at location "
                + getMem(*this));
    //##########################################################################
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Constructor with name and UGrid
  //----------------------------------------------------------------------------
  template<typename T>
  VectorField<T>::VectorField(std::string name,
                              std::shared_ptr<UGrid<T>> ugrid)
  : _name(name), _ugrid(ugrid)
  {
    _N = _ugrid->getN();
    _dim = _ugrid->getDim();
    _approx = std::make_shared<Interpolator<T>>();
    //##########################################################################
    _log = std::make_shared<Log>();
    _log->init("ET:VectorField:" + _name, ".logs/vectorfield_default.txt");
    _log->TRACE("Vector Field '" + _name + "' created at location "
                + getMem(*this));
    //##########################################################################
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Constructor with UGrid and field
  //----------------------------------------------------------------------------
  template<typename T>
  VectorField<T>::VectorField(std::shared_ptr<UGrid<T>> ugrid,
                              std::vector<std::vector<T>> field)
  : _ugrid(ugrid), _field(field), _name("default")
  {
    _N = ugrid->getN();
    _dim = _ugrid->getDim();
    _approx = std::make_shared<Interpolator<T>>();
    //##########################################################################
    _log = std::make_shared<Log>();
    _log->init("ET:VectorField:default", ".logs/vectorfield_default.txt");
    _log->TRACE("Vector Field 'default' created at location "
                + getMem(*this));
    //##########################################################################
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Constructor with name, UGrid and field
  //----------------------------------------------------------------------------
  template<typename T>
  VectorField<T>::VectorField(std::string name,
                              std::shared_ptr<UGrid<T>> ugrid,
                              std::vector<std::vector<T>> field)
  : _name(name), _ugrid(ugrid), _field(field)
  {
    _N = ugrid->getN();
    _dim = _ugrid->getDim();
    _approx = std::make_shared<Interpolator<T>>();
    //##########################################################################
    _log = std::make_shared<Log>();
    _log->init("ET:VectorField:" + _name, ".logs/vectorfield_default.txt");
    _log->TRACE("Vector Field '" + _name + "' created at location "
                + getMem(*this));
    //##########################################################################
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //	Various constructors taking in shared logger
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //	Constructor with shared logger
  //----------------------------------------------------------------------------
  template<typename T>
  VectorField<T>::VectorField(std::shared_ptr<Log> log)
  : _name("default"), _dim(0), _N(0)
  {
    //##########################################################################
    _log = log;
    _log->TRACE("Vector Field 'default' created at location "
                + getMem(*this));
    _log->INFO("Logger passed to Vector Field 'default'");
    //##########################################################################
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //	Constructor with UGrid and shared logger
  //----------------------------------------------------------------------------
  template<typename T>
  VectorField<T>::VectorField(std::shared_ptr<UGrid<T>> ugrid,
                              std::shared_ptr<Log> log)
  : _ugrid(ugrid), _name("default")
  {
    _N = _ugrid->getN();
    _dim = _ugrid->getDim();
    _approx = std::make_shared<Interpolator<T>>();
    //##########################################################################
    _log = log;
    _log->TRACE("Vector Field 'default' created at location "
                + getMem(*this));
    _log->INFO("Logger passed to Vector Field 'default'");
    //##########################################################################
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Constructor with name and UGrid
  //----------------------------------------------------------------------------
  template<typename T>
  VectorField<T>::VectorField(std::string name,
                              std::shared_ptr<UGrid<T>> ugrid,
                              std::shared_ptr<Log> log)
  : _name(name), _ugrid(ugrid)
  {
    _N = _ugrid->getN();
    _dim = _ugrid->getDim();
    _approx = std::make_shared<Interpolator<T>>();
    //##########################################################################
    _log = log;
    _log->TRACE("Vector Field '" + _name + "' created at location "
                + getMem(*this));
    _log->INFO("Logger passed to Vector Field '" + _name + "'");
    //##########################################################################
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Constructor with UGrid and field
  //----------------------------------------------------------------------------
  template<typename T>
  VectorField<T>::VectorField(std::shared_ptr<UGrid<T>> ugrid,
                              std::vector<std::vector<T>> field,
                              std::shared_ptr<Log> log)
  : _ugrid(ugrid), _field(field), _name("default")
  {
    _N = ugrid->getN();
    _dim = _ugrid->getDim();
    _approx = std::make_shared<Interpolator<T>>();
    //##########################################################################
    _log = log;
    _log->TRACE("Vector Field 'default' created at location "
                + getMem(*this));
    _log->INFO("Logger passed to Vector Field 'default'");
    //##########################################################################
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Constructor with name, UGrid and field
  //----------------------------------------------------------------------------
  template<typename T>
  VectorField<T>::VectorField(std::string name,
                              std::shared_ptr<UGrid<T>> ugrid,
                              std::vector<std::vector<T>> field,
                              std::shared_ptr<Log> log)
  : _name(name), _ugrid(ugrid), _field(field)
  {
    _N = ugrid->getN();
    _dim = _ugrid->getDim();
    _approx = std::make_shared<Interpolator<T>>();
    //##########################################################################
    _log = log;
    _log->TRACE("Vector Field '" + _name + "' created at location "
                + getMem(*this));
    _log->INFO("Logger passed to Vector Field '" + _name + "'");
    //##########################################################################
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Getters and Setters
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<std::vector<T>> VectorField<T>::getField() const
  {
    return _field;
  }
  template<typename T>
  std::vector<std::vector<T>>* VectorField<T>::accessField()
  {
    return &_field;
  }
  template<typename T>
  T* VectorField<T>::data()
  {
    return _field[0].data();
  }
  template<typename T>
  std::string VectorField<T>::getName() const
  {
    return _name;
  }
  template<typename T>
  uint64_t VectorField<T>::getN() const
  {
    return _N;
  }
  template<typename T>
  uint32_t VectorField<T>::getDim() const
  {
    return _dim;
  }
  template<typename T>
  std::shared_ptr<Interpolator<T>> VectorField<T>::getInterpolator() const
  {
    return _approx;
  }
  template<typename T>
  int VectorField<T>::getFlag() const
  {
    return _flag;
  }
  template<typename T>
  std::string VectorField<T>::getInfo() const
  {
    return _info;
  }
  template<typename T>
  std::shared_ptr<Log> VectorField<T>::getLogger()
  {
    return _log;
  }
  template<typename T>
  void VectorField<T>::setUGrid(std::shared_ptr<UGrid<T>> ugrid)
  {
    _ugrid = ugrid;
  }
  template<typename T>
  void VectorField<T>::setField(std::vector<std::vector<T>> field)
  {
    _field = field;
  }
  template<typename T>
  void VectorField<T>::setName(std::string name)
  {
    _name = name;
  }
  template<typename T>
  void VectorField<T>::setInterpolatorType(std::string type)
  {
    _approx->setInterpolatorType(type);
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Operator overloads
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<T>& VectorField<T>::operator()(const uint32_t& i)
  {
    if (i >= _N )
    {
      //########################################################################
      _log->ERROR("VectorField " + _name
                  + ": Attempted to access _field array of size "
                  + std::to_string(_N) + " with index "
                  + std::to_string(i));
      //########################################################################
      if(_field.size() > 0)
      {
        //######################################################################
        _log->INFO("VectorField "+ _name +": Returning the element at index 0");
        //######################################################################
        return _field[0];
      }
      else
      {
        //######################################################################
        _log->INFO("VectorField " + _name + ": Terminating program");
        //######################################################################
        exit(0);
      }
    }
    return _field[i];
  }
  template<typename T>
  const std::vector<T>& VectorField<T>::operator()(const uint32_t& i) const
  {
    if (i >= _N )
    {
      //########################################################################
      _log->ERROR("VectorField " + _name
                  + ": Attempted to access _field array of size "
                  + std::to_string(_N) + " with index "
                  + std::to_string(i));
      //########################################################################
      if(_field.size() > 0)
      {
        //######################################################################
        _log->INFO("VectorField "+ _name +": Returning the element at index 0");
        //######################################################################
        return _field[0];
      }
      else
      {
        //######################################################################
        _log->INFO("VectorField " + _name + ": Terminating program");
        //######################################################################
        exit(0);
      }
    }
    return _field[i];
  }
  template<typename T>
  T& VectorField<T>::operator()(const uint32_t& i, const uint32_t& j)
  {
    return _field[i][j];
  }
  template<typename T>
  const T& VectorField<T>::operator()(const uint32_t& i, const uint32_t& j) const
  {
    return _field[i][j];
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Derivative operations
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Derivative of a vector field along a particular direction
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<std::vector<T>>
  VectorField<T>::derivative(uint32_t dir, uint32_t n)
  {
    return _approx->vectorDerivative(_ugrid, (*this), dir, n);
  }
  //----------------------------------------------------------------------------
}
