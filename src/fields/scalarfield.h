//------------------------------------------------------------------------------
//  scalarfield.h
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
  //  Scalar fields class.
  //----------------------------------------------------------------------------
  template<typename T>
  class ScalarField: public std::enable_shared_from_this<ScalarField<T>>
  {
  public:
    //--------------------------------------------------------------------------
    //  Constructors
    //--------------------------------------------------------------------------
    ScalarField();
    virtual ~ScalarField();
    ScalarField(std::shared_ptr<UGrid<T>> ugrid);
    ScalarField(std::string, std::shared_ptr<UGrid<T>> ugrid);
    ScalarField(std::shared_ptr<UGrid<T>> ugrid, std::vector<T> field);
    ScalarField(std::string, std::shared_ptr<UGrid<T>> ugrid,
                std::vector<T> field);
    //--------------------------------------------------------------------------
    //  Constructors with shared loggers
    //--------------------------------------------------------------------------
    ScalarField(std::shared_ptr<Log> log);
    ScalarField(std::shared_ptr<UGrid<T>> ugrid, std::shared_ptr<Log> log);
    ScalarField(std::string, std::shared_ptr<UGrid<T>> ugrid,
                std::shared_ptr<Log> log);
    ScalarField(std::shared_ptr<UGrid<T>> ugrid, std::vector<T> field,
                std::shared_ptr<Log> log);
    ScalarField(std::string, std::shared_ptr<UGrid<T>> ugrid,
                std::vector<T> field, std::shared_ptr<Log> log);
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  Getters
    //--------------------------------------------------------------------------
    std::vector<T> getField() const;
    std::vector<T>* accessField();
    T* data();
    std::string getName() const;
    uint64_t getN() const;
    uint32_t getDim() const;
    std::shared_ptr<Approximator<T>> getApproximator() const;
    std::shared_ptr<UGrid<T>> getUGrid() const;
    std::shared_ptr<Log> getLogger();
    int getFlag() const;
    std::string getInfo() const;
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  Setters
    //--------------------------------------------------------------------------
    void setUGrid(std::shared_ptr<UGrid<T>> ugrid);
    void setField(std::vector<T> field);
    void setName(std::string name);
    void setApproxType(std::string type);
    void setFlag(int flag);
    void setInfo(std::string info);
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  Operator overloads
    //--------------------------------------------------------------------------
    ScalarField<T> operator+(const ScalarField<T>& scalar) const;
    ScalarField<T> operator-(const ScalarField<T>& scalar) const;
    ScalarField<T> operator*(const ScalarField<T>& scalar) const;
    ScalarField<T> operator/(const ScalarField<T>& scalar) const;
    ScalarField<T>& operator+=(const ScalarField<T>& scalar);
    ScalarField<T>& operator-=(const ScalarField<T>& scalar);
    ScalarField<T>& operator*=(const ScalarField<T>& scalar);
    ScalarField<T>& operator/=(const ScalarField<T>& scalar);
    T& operator()(const uint32_t& i);
    const T& operator()(const uint32_t& i) const;
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  Methods for calculating derivatives
    //--------------------------------------------------------------------------
    Matrix<T> constructTaylorMatrix();
    std::vector<std::vector<T>> gradient();
    std::vector<std::vector<T>> derivative(uint32_t n);
    std::vector<T> derivative(uint32_t dir, uint32_t n);
    std::vector<T> derivative(std::vector<uint32_t> deriv);
    std::vector<T> derivativePoint(uint64_t index, uint32_t n);
    T derivativePoint(uint64_t index, uint32_t dir, uint32_t n);
    T derivativePoint(uint64_t index, std::vector<uint32_t> deriv);
    std::vector<T> laplacian();
    T laplacian(uint64_t index);
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  Various functions
    //--------------------------------------------------------------------------
    std::string summary();
    //--------------------------------------------------------------------------

  private:
    //--------------------------------------------------------------------------
    //  Basic attributes
    //--------------------------------------------------------------------------
    std::string _name;
    uint32_t _dim;
    uint64_t _N;
    std::vector<T> _field;
    //--------------------------------------------------------------------------
    //  Shared objects
    //--------------------------------------------------------------------------
    std::shared_ptr<UGrid<T>> _ugrid;
    std::shared_ptr<Approximator<T>> _approx;
    std::shared_ptr<Log> _log;

    int _flag;
    std::string _info;
    //--------------------------------------------------------------------------

  protected:
    //--------------------------------------------------------------------------
    // Inherited functions which much be overwritten
    //--------------------------------------------------------------------------
    std::vector<std::vector<T>> diffEQ();
    std::vector<T> diffEQ(uint64_t index);
    //--------------------------------------------------------------------------

  };
  //----------------------------------------------------------------------------

  //  Explicit instantiation of double
  template class ScalarField<double>;

  //----------------------------------------------------------------------------
  //  Example scalar fields
  //----------------------------------------------------------------------------
  template<typename T>
  class KleinGordon1D : public ScalarField<T>
  {
  public:
    KleinGordon1D();
    KleinGordon1D(std::shared_ptr<UGrid<T>> ugrid);
    KleinGordon1D(std::shared_ptr<UGrid<T>> ugrid, T mass);
    //  Getters
    T getMass();
    //  Setters
    void setMass(T mass);
  private:
    T _mass;
  protected:
    std::vector<std::vector<T>> diffEQ();
    std::vector<T> diffEQ(uint64_t index);
  };
  //----------------------------------------------------------------------------

  template class KleinGordon1D<double>;

}
