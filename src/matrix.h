//  Matrix class for Entropic Trajectories and EDMC
#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <stdio.h>
#include <complex>
#include <lapacke.h>
#include <cblas.h>
#include "vector.h"

namespace ET
{
  template<typename T> class Matrix
  {
  public:
    Matrix();
    ~Matrix();
    Matrix(const Matrix<T>& m);
    Matrix(unsigned int n);
    Matrix(std::string name, unsigned int n);
    Matrix(unsigned int n, unsigned int m);
    Matrix(std::string name, unsigned int n, unsigned int m);
    Matrix(unsigned int n, unsigned int m, const T& init);
    Matrix(std::string name, unsigned int n, unsigned int m, const T& init);

    //  Constructors passing elements
    Matrix(std::vector<std::vector<T> > mat);
    Matrix(std::string name, std::vector<std::vector<T> > mat);
    Matrix(unsigned int n, std::vector<T> flat);
    Matrix(std::string name, unsigned int n, std::vector<T> flat);
    Matrix(unsigned int n, unsigned int m, std::vector<T> flat);
    Matrix(std::string name, unsigned int n,
      unsigned int m, std::vector<T> flat);
    //  Getters
    unsigned int get_rows() const;
    unsigned int get_cols() const;
    std::string get_name() const;

    //  Setters
    void set_name(std::string name);

    //  Operator overloads
    Matrix<T>& operator=(const Matrix<T>& m);
    bool operator==(const Matrix<T>& m);
    //  Matrix algebra
    Matrix<T> operator+(const Matrix<T>& m);
    Matrix<T>& operator+=(const Matrix<T>& m);
    Matrix<T> operator-(const Matrix<T>& m);
    Matrix<T>& operator-=(const Matrix<T>& m);
    Matrix<T> operator*(const Matrix<T>& m);
    Matrix<T>& operator*=(const Matrix<T>& m);
    //  Scalar algebra
    Matrix<T> operator+(const T& s);
    Matrix<T> operator-(const T& s);
    Matrix<T> operator*(const T& s);
    Matrix<T> operator/(const T& s);
    Matrix<T>& operator+=(const T& s);
    Matrix<T>& operator-=(const T& s);
    Matrix<T>& operator*=(const T& s);
    Matrix<T>& operator/=(const T& s);
    //  Multiplying a vector
    Vector<T> operator*(const Vector<T>& v);
    //  Access operators
    T& operator()(const unsigned int& i, const unsigned int& j);
    const T& operator()(const unsigned int& i, const unsigned int& j) const;

    //  Various methods
    void print();
    Matrix<T> transpose();

  private:
    //  _n is the number of rows, _m is the number of columns
    unsigned int _n, _m;
    std::vector<std::vector<T> > _mat;
    //  possible name for the matrix
    std::string _name;
  };

  //  Various matrices
  template<typename T>
  Matrix<T> identity(unsigned int n);
  template<typename T>
  Matrix<T> zeroes(unsigned int n);
  template<typename T>
  Matrix<T> zeroes(unsigned int n, unsigned int m);
  template<typename T>
  Matrix<T> ones(unsigned int n);
  template<typename T>
  Matrix<T> ones(unsigned int n, unsigned int m);

}
