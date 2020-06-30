//------------------------------------------------------------------------------
//  scalarfield.h
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
#include <memory>

#include "ugrid.h"

namespace ET
{
  template<typename T> class Approximator;
}
#include "approximator.h"
namespace ET
{
  template<typename T>
  class ScalarField
  {
  public:
    ScalarField();
    ~ScalarField();
    ScalarField(std::shared_ptr<UGrid<T>> micro);
    ScalarField(std::string, std::shared_ptr<UGrid<T>> micro);
    ScalarField(std::shared_ptr<UGrid<T>> micro, std::vector<T> field);
    ScalarField(std::string, std::shared_ptr<UGrid<T>> micro, std::vector<T> field);

    //  Getters
    //  get const reference to field
    std::vector<T> getField() const;
    //  get access to field
    std::vector<T>* accessField();
    //  get access to beginning of field
    T* data();
    std::string getName() const;
    uint64_t getN() const;
    Approximator<T>* getApproximator();
    int getFlag();
    std::string getInfo();

    //  Setters
    void setUGrid(std::shared_ptr<UGrid<T>> micro);
    void setField(std::vector<T> field);
    void setName(std::string name);
    void setApproxType(std::string type);
    void setFlag(int flag);
    void setInfo(std::string info);

    //  operator overloads
    T& operator()(const uint32_t& i);
    const T& operator()(const uint32_t& i) const;

    //  Methods for calculating derivatives
    Matrix<T> constructTaylorMatrix();


  private:
    //  number of points
    uint64_t _N;
    //  pointer to associated microstates
    std::shared_ptr<UGrid<T>> _micro;
    //  vector for field values
    std::vector<T> _field;
    //  pointer to associated approximator
    Approximator<T>* _approx;
    //  name of the field
    std::string _name;
    //  conatiner for message status
    int _flag;
    //  container for messages
    std::string _info;

  };

  template class ScalarField<double>;
}
