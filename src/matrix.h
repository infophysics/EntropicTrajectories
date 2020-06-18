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
#include <stdint.h>

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
    Matrix(uint64_t n);
    Matrix(std::string name, uint64_t n);
    Matrix(uint64_t n, uint64_t m);
    Matrix(std::string name, uint64_t n, uint64_t m);
    Matrix(uint64_t n, uint64_t m, const T& init);
    Matrix(std::string name, uint64_t n, uint64_t m, const T& init);

    //  Constructors passing elements
    Matrix(uint64_t n, std::vector<T> flat);
    Matrix(std::string name, uint64_t n, std::vector<T> flat);
    Matrix(uint64_t n, uint64_t m, std::vector<T> flat);
    Matrix(std::string name, uint64_t n,
      uint64_t m, std::vector<T> flat);
    Matrix(std::vector<std::vector<T> > array);
    Matrix(std::string name, std::vector<std::vector<T> > array);

    //  Getters
    uint64_t getNumRows() const;
    uint64_t getNumCols() const;
    std::string getName() const;
    std::vector<T> getArray() const;
    std::vector<T> getRow(uint64_t i);
    std::vector<T> getCol(uint64_t i);

    //  Setters
    void setName(std::string name);
    void setRow(uint64_t i, std::vector<T> row);
    void setCol(uint64_t i, std::vector<T> col);
    void setArray(uint64_t m, std::vector<T> mat);
    void setArray(std::vector<std::vector<T> > mat);

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
      uint64_t n = matrix.getNumRows();
      uint64_t m = matrix.getNumCols();
      std::string name = "(" + std::to_string(s) + "I + "  + matrix.getName() + ")";
      Matrix<T> l(name,n,m,0.0);
      for (uint64_t i = 0; i < n*m; i++) {
          l(i) = matrix(i) + s;
      }
      return l;
    }
    friend Matrix<T> operator-(T s, const Matrix<T>& matrix)
    {
      uint64_t n = matrix.getNumRows();
      uint64_t m = matrix.getNumCols();
      std::string name = "(" + std::to_string(s) + "I - "  + matrix.getName() + ")";
      Matrix<T> l(name,n,m,0.0);
      for (uint64_t i = 0; i < n*m; i++) {
          l(i) = s - matrix(i);
      }
      return l;
    }
    friend Matrix<T> operator*(T s, const Matrix<T>& matrix)
    {
      uint64_t n = matrix.getNumRows();
      uint64_t m = matrix.getNumCols();
      std::string name = "(" + std::to_string(s) + " * "  + matrix.getName() + ")";
      Matrix<T> l(name,n,m,0.0);
      for (uint64_t i = 0; i < n*m; i++) {
          l(i) = matrix(i) * s;
      }
      return l;
    }
    friend Matrix<T> operator/(T s, const Matrix<T>& matrix)
    {
      uint64_t n = matrix.getNumRows();
      uint64_t m = matrix.getNumCols();
      std::string name = "(" + std::to_string(s) + " / "  + matrix.getName() + ")";
      Matrix<T> l(name,n,m,0.0);
      std::vector<T> mat(n*m);
      for (uint64_t i = 0; i < n*m; i++)
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
      l.setArray(m, mat);
      return l;
    }
    //  Multiplying a vector
    Vector<T> operator*(const Vector<T>& v);
    //  Access operators
    T& operator()(const uint64_t& i, const uint64_t& j);
    const T& operator()(const uint64_t& i, const uint64_t& j) const;
    //  Flattened access
    T& operator()(const uint64_t& i);
    const T& operator()(const uint64_t& i) const;

    //  Linear algebra tools
    Matrix<T> permutationMatrix(int& n, std::vector<int>& pivot);
    Matrix<T> inverse();
    Matrix<T> LU();
    void QR();
    void SVD();
    std::vector<T> singularValues(int info);

    //  Various methods
    void print();
    std::string summary();
    Matrix<T> transpose();
    void transpose(bool inplace = true);

  private:
    //  _n is the number of rows, _m is the number of columns
    uint64_t _n, _m;
    std::vector<T> _array;
    //  possible name for the matrix
    std::string _name;
  };

  //  Various matrices
  template<typename T>
  Matrix<T> identity(uint64_t n);
  template<typename T>
  Matrix<T> zeroes(uint64_t n);
  template<typename T>
  Matrix<T> zeroes(uint64_t n, uint64_t m);
  template<typename T>
  Matrix<T> ones(uint64_t n);
  template<typename T>
  Matrix<T> ones(uint64_t n, uint64_t m);

  template<typename T>
  std::ostream& operator<<(std::ostream& os, const Matrix<T>& matrix);

  template class Matrix<double>;
}
