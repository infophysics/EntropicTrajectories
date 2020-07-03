//------------------------------------------------------------------------------
//  approximator.h
//  The Entropic Trajectories Framework
//  -----------------------------------
//  Copyright (C) [2020] by [N. Carrara, F. Costa, P. Pessoa]
//  [ncarrara@albany.edu,felipecosta.physics@gmail.com,
//    pedroh.pessoa100@gmail.com]
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
#include "params.h"
#include "utils.h"
#include "matrix.h"

//------------------------------------------------------------------------------
//  Forward declaration of ScalarField
//------------------------------------------------------------------------------
namespace ET
{
  template<typename T> class ScalarField;
}
#include "scalarfield.h"
namespace ET
{
  //----------------------------------------------------------------------------
  //  Approximator class.
  //----------------------------------------------------------------------------
  template<typename T>
  class Approximator
  {
  public:
    //--------------------------------------------------------------------------
    //  Constructors
    //--------------------------------------------------------------------------
    Approximator();
    ~Approximator();
    Approximator(int type);
    Approximator(std::string type);
    //--------------------------------------------------------------------------
    //  Constructors with shared loggers
    //--------------------------------------------------------------------------
    Approximator(std::shared_ptr<Log> log);
    Approximator(int type, std::shared_ptr<Log> log);
    Approximator(std::string type, std::shared_ptr<Log> log);
    //--------------------------------------------------------------------------
    //  Getters
    //--------------------------------------------------------------------------
    int getApproxType() const;
    ApproxParams getApproxParams() const;
    int getLSDriver() const;
    int getFlag() const;
    std::string getInfo() const;
    std::shared_ptr<Log> getLogger();
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  Setters
    //--------------------------------------------------------------------------
    void setApproxType(std::string type);
    void setApproxParams(ApproxParams params);
    void setLSDriver(std::string type);
    //  Set parameters
    void set_k(uint64_t k);
    void set_n(uint64_t n);
    void setFlag(int flag);
    void setInfo(std::string info);
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  Gradient functions
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  Scalar fields
    //--------------------------------------------------------------------------
    std::vector<T>
    scalarGradientPoint(const std::shared_ptr<UGrid<T>> ugrid,
                        const std::shared_ptr<ScalarField<T>> field,
                        uint64_t index);
    std::vector<T>
    scalarGradientMLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
                           const std::shared_ptr<ScalarField<T>> field,
                           uint64_t index);
    std::vector<std::vector<T>>
    scalarGradient(const std::shared_ptr<UGrid<T>> ugrid,
                   const std::shared_ptr<ScalarField<T>> field);
    std::vector<std::vector<T>>
    scalarGradientMLS(const std::shared_ptr<UGrid<T>> ugrid,
                      const std::shared_ptr<ScalarField<T>> field);
    //--------------------------------------------------------------------------
    //  Passing field as a const refernce
    //--------------------------------------------------------------------------
    std::vector<T>
    scalarGradientPoint(const std::shared_ptr<UGrid<T>> ugrid,
                        const ScalarField<T>& field,
                        uint64_t index);
    std::vector<T>
    scalarGradientMLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
                           const ScalarField<T>& field,
                           uint64_t index);
    std::vector<std::vector<T>>
    scalarGradient(const std::shared_ptr<UGrid<T>> ugrid,
                   const ScalarField<T>& field);
    std::vector<std::vector<T>>
    scalarGradientMLS(const std::shared_ptr<UGrid<T>> ugrid,
                      const ScalarField<T>& field);
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  nth-derivatives of scalar field
    //--------------------------------------------------------------------------
    std::vector<T> scalarDerivative(const std::shared_ptr<UGrid<T>> ugrid,
                                    const std::shared_ptr<ScalarField<T>> field,
                                    uint32_t dir, uint32_t n);
    //--------------------------------------------------------------------------
    //  Passing field as a const refernce
    //--------------------------------------------------------------------------
    std::vector<T> scalarDerivative(const std::shared_ptr<UGrid<T>> ugrid,
                                    const ScalarField<T>& field,
                                    uint32_t dir, uint32_t n);
    //--------------------------------------------------------------------------


    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  MLS functions
    //--------------------------------------------------------------------------
    Matrix<T> constructTaylorMatrix(const std::shared_ptr<UGrid<T>> ugrid,
                                    const std::vector<uint64_t> neighbors,
                                    uint64_t index, uint64_t order);
    Matrix<T> constructTaylorMatrix(const std::shared_ptr<UGrid<T>> ugrid,
                                    const std::vector<uint64_t> neighbors,
                                    uint64_t index, Monomial& mono);
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  various functions
    //--------------------------------------------------------------------------
    std::string summary();
    //--------------------------------------------------------------------------

  private:
    //--------------------------------------------------------------------------
    //  Basic attributes
    //--------------------------------------------------------------------------
    std::string _name;
    enum ApproxType _type;
    struct ApproxParams _params;
    enum LSDriver _lsdriver;
    //--------------------------------------------------------------------------
    //  Shared objects
    //--------------------------------------------------------------------------
    std::shared_ptr<Log> _log;

    //  conatiner for message status
    int _flag;
    //  container for messages
    std::string _info;
  };


  template class Approximator<double>;
}
