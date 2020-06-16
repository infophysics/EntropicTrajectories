//  Matrix class for Entropic Trajectories and EDMC
#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <ostream>
#include <stdio.h>
#include <complex>
#include <lapacke.h>
#include <cblas.h>
#include "vector.h"
#include "utils.h"

namespace ET
{
  //  This matrix class acts as a general container for n x m matrices
  //  which are row ordered, i.e. m(i,j) picks the element of the ith row
  //  and the jth column.  It wraps several methods from BLAS and LAPACK.
  //
  template<typename T>
  class Matrix
  {
  public:
    Matrix();
    ~Matrix();
    Matrix(const Matrix<T>& matrix);
    Matrix(unsigned int n);
    Matrix(std::string name, unsigned int n);
    Matrix(unsigned int n, unsigned int m);
    Matrix(std::string name, unsigned int n, unsigned int m);
    Matrix(unsigned int n, unsigned int m, const T& init);
    Matrix(std::string name, unsigned int n, unsigned int m, const T& init);

    //  Constructors passing elements
    Matrix(unsigned int n, std::vector<T> flat);
    Matrix(std::string name, unsigned int n, std::vector<T> flat);
    Matrix(unsigned int n, unsigned int m, std::vector<T> flat);
    Matrix(std::string name, unsigned int n,
      unsigned int m, std::vector<T> flat);
    Matrix(std::vector<std::vector<T> > array);
    Matrix(std::string name, std::vector<std::vector<T> > array);
    //  Getters
    unsigned int get_rows() const;
    unsigned int get_cols() const;
    std::string get_name() const;
    std::vector<T> get_mat() const;
    std::vector<T> get_row(unsigned int i);
    std::vector<T> get_col(unsigned int i);


    //  Setters
    void set_name(std::string name);
    void set_row(unsigned int i, std::vector<T> row);
    void set_col(unsigned int i, std::vector<T> col);
    void set_mat(unsigned int m, std::vector<T> mat);
    void set_mat(std::vector<std::vector<T> > mat);

    //  Operator overloads
    Matrix<T>& operator=(const Matrix<T>& matrix);
    bool operator==(const Matrix<T>& matrix) const;
    bool operator!=(const Matrix<T>& matrix) const;
    Matrix<T> operator-() const;
    //  Matrix algebra
    Matrix<T> operator+(const Matrix<T>& matrix) const;
    Matrix<T>& operator+=(const Matrix<T>& matrix);
    Matrix<T> operator-(const Matrix<T>& matrix) const;
    Matrix<T>& operator-=(const Matrix<T>& matrix);
    //  Using CBLAS for multiplication
    Matrix<T> operator*(const Matrix<T>& matrix) const;
    Matrix<T> brute_mul(const Matrix<T>& matrix) const;
    Matrix<T>& operator*=(const Matrix<T>& matrix);
    //  Scalar algebra
    Matrix<T> operator+(const T& s) const;
    Matrix<T> operator-(const T& s) const;
    Matrix<T> operator*(const T& s) const;
    Matrix<T> operator/(const T& s) const;
    Matrix<T>& operator+=(const T& s);
    Matrix<T>& operator-=(const T& s);
    Matrix<T>& operator*=(const T& s);
    Matrix<T>& operator/=(const T& s);
    //  Overloads of scalar operations from the left.
    //  since we are trying to friend a template argument,
    //  the friend method must be defined within the class block.
    friend Matrix<T> operator+(T s, const Matrix<T>& matrix)
    {
      unsigned int n = matrix.get_rows();
      unsigned int m = matrix.get_cols();
      Matrix<T> l(n,m,0.0);
      for (unsigned int i = 0; i < n*m; i++) {
          l(i) = matrix(i) + s;
      }
      return l;
    }
    friend Matrix<T> operator-(T s, const Matrix<T>& matrix)
    {
      unsigned int n = matrix.get_rows();
      unsigned int m = matrix.get_cols();
      Matrix<T> l(n,m,0.0);
      for (unsigned int i = 0; i < n*m; i++) {
          l(i) = s - matrix(i);
      }
      return l;
    }
    friend Matrix<T> operator*(T s, const Matrix<T>& matrix)
    {
      unsigned int n = matrix.get_rows();
      unsigned int m = matrix.get_cols();
      Matrix<T> l(n,m,0.0);
      for (unsigned int i = 0; i < n*m; i++) {
          l(i) = matrix(i) * s;
      }
      return l;
    }
    friend Matrix<T> operator/(T s, const Matrix<T>& matrix)
    {
      unsigned int n = matrix.get_rows();
      unsigned int m = matrix.get_cols();
      Matrix<T> l(n,m,0.0);
      std::vector<T> mat(n*m);
      for (unsigned int i = 0; i < n*m; i++)
      {
        if (matrix(i) == 0)
        {
          return l;
        }
        else
        {
          mat[i] = s / matrix(i);
        }
      }
      l.set_mat(m, mat);
      return l;
    }
    //  Multiplying a vector
    Vector<T> operator*(const Vector<T>& v);
    //  Access operators
    T& operator()(const unsigned int& i, const unsigned int& j);
    const T& operator()(const unsigned int& i, const unsigned int& j) const;
    //  Flattened access
    T& operator()(const unsigned int& i);
    const T& operator()(const unsigned int& i) const;

    //  Various methods
    void print();
    std::string summary();
    Matrix<T> transpose();
    void transpose(bool inplace = true);

  private:
    //  _n is the number of rows, _m is the number of columns
    unsigned int _n, _m;
    std::vector<T> _mat;
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

  template<typename T>
  std::ostream& operator<<(std::ostream& os, const Matrix<T>& matrix);

  template class Matrix<double>;
}
