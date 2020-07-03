//------------------------------------------------------------------------------
//  vectorfield.h
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

#include "ugrid.h"
#include "log.h"

//------------------------------------------------------------------------------
//  Forward declaration of Approximator
//------------------------------------------------------------------------------
namespace ET
{
  template<typename T> class Approximator;
}
#include "approximator.h"
namespace ET
{
  //----------------------------------------------------------------------------
  //  Vector fields class.
  //----------------------------------------------------------------------------
  template<typename T>
  class VectorField: public std::enable_shared_from_this<VectorField<T>>
  {
  public:
    //--------------------------------------------------------------------------
    //  Constructors
    //--------------------------------------------------------------------------
    VectorField();
    ~VectorField();
    VectorField(std::shared_ptr<UGrid<T>> ugrid);
    VectorField(std::string, std::shared_ptr<UGrid<T>> ugrid);
    VectorField(std::shared_ptr<UGrid<T>> ugrid,
                std::vector<std::vector<T>> field);
    VectorField(std::string, std::shared_ptr<UGrid<T>> ugrid,
                std::vector<T> field);
    //--------------------------------------------------------------------------
    //  Constructors with shared loggers
    //--------------------------------------------------------------------------
    VectorField(std::shared_ptr<Log> log);
    VectorField(std::shared_ptr<UGrid<T>> ugrid, std::shared_ptr<Log> log);
    VectorField(std::string, std::shared_ptr<UGrid<T>> ugrid,
                std::shared_ptr<Log> log);
    VectorField(std::shared_ptr<UGrid<T>> ugrid,
                std::vector<std::vector<T>> field,
                std::shared_ptr<Log> log);
    VectorField(std::string, std::shared_ptr<UGrid<T>> ugrid,
                std::vector<T> field, std::shared_ptr<Log> log);
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  Getters
    //--------------------------------------------------------------------------
    std::vector<std::vector<T>> getField() const;
    std::vector<std::vector<T>>* accessField();
    T* data();
    std::string getName() const;
    uint64_t getN() const;
    uint32_t getDim() const;
    std::shared_ptr<Approximator<T>> getApproximator() const;
    int getFlag() const;
    std::string getInfo() const;
    std::shared_ptr<Log> getLogger();
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  Setters
    //--------------------------------------------------------------------------
    void setUGrid(std::shared_ptr<UGrid<T>> ugrid);
    void setField(std::vector<std::vector<T>> field);
    void setName(std::string name);
    void setApproxType(std::string type);
    void setFlag(int flag);
    void setInfo(std::string info);
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  Operator overloads
    //--------------------------------------------------------------------------
    T& operator()(const uint32_t& i);
    const T& operator()(const uint32_t& i) const;
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  Methods for calculating derivatives
    //--------------------------------------------------------------------------
    Matrix<T> constructTaylorMatrix();
    std::vector<std::vector<T>> gradient();
    std::vector<T> derivative(uint32_t dir, uint32_t n);
    T divergence();
    //--------------------------------------------------------------------------

  private:
    //--------------------------------------------------------------------------
    //  Basic attributes
    //--------------------------------------------------------------------------
    std::string _name;
    uint32_t _dim;
    uint64_t _N;
    std::vector<std::vector<T>> _field;
    //--------------------------------------------------------------------------
    //  Shared objects
    //--------------------------------------------------------------------------
    std::shared_ptr<UGrid<T>> _ugrid;
    std::shared_ptr<Approximator<T>> _approx;
    std::shared_ptr<Log> _log;

    int _flag;
    std::string _info;
    //--------------------------------------------------------------------------

  };
  //----------------------------------------------------------------------------

  //  Explicit instantiation of double
  template class VectorField<double>;
}
