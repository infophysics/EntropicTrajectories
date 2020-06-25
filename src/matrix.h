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
#include <complex>

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
    Matrix(uint32_t m);
    Matrix(std::string name, uint32_t m);
    Matrix(uint32_t m, uint32_t n);
    Matrix(std::string name, uint32_t m, uint32_t n);
    Matrix(uint32_t m, uint32_t n, const T& init);
    Matrix(std::string name, uint32_t m, uint32_t n, const T& init);

    //  Constructors passing elements
    Matrix(uint32_t m, std::vector<T> flat);
    Matrix(std::string name, uint32_t m, std::vector<T> flat);
    Matrix(uint32_t m, uint32_t n, std::vector<T> flat);
    Matrix(std::string name, uint32_t m,
      uint32_t n, std::vector<T> flat);
    Matrix(std::string name, uint32_t m, uint32_t n, T* array);
    Matrix(std::vector<std::vector<T> > array);
    Matrix(std::string name, std::vector<std::vector<T> > array);

    //  Getters
    uint32_t getNumRows() const;
    uint32_t getNumCols() const;
    std::string getName() const;
    //  get a const reference to array
    std::vector<T> getArray() const;
    //  get access to array
    std::vector<T> accessArray();
    std::vector<T> getRow(uint32_t i);
    std::vector<T> getCol(uint32_t i);
    float *data();

    //  Setters
    void setName(std::string name);
    void setRow(uint32_t i, std::vector<T> row);
    void setCol(uint32_t i, std::vector<T> col);
    void setArray(uint32_t m, std::vector<T> mat);
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
      uint32_t m = matrix.getNumRows();
      uint32_t n = matrix.getNumCols();
      std::string name = "(" + std::to_string(s) + "I + "  + matrix.getName() + ")";
      Matrix<T> l(name,m,n,0.0);
      for (uint32_t i = 0; i < m*n; i++) {
          l(i) = matrix(i) + s;
      }
      return l;
    }
    friend Matrix<T> operator-(T s, const Matrix<T>& matrix)
    {
      uint32_t m = matrix.getNumRows();
      uint32_t n = matrix.getNumCols();
      std::string name = "(" + std::to_string(s) + "I - "  + matrix.getName() + ")";
      Matrix<T> l(name,m,n,0.0);
      for (uint32_t i = 0; i < m*n; i++) {
          l(i) = s - matrix(i);
      }
      return l;
    }
    friend Matrix<T> operator*(T s, const Matrix<T>& matrix)
    {
      uint32_t m = matrix.getNumRows();
      uint32_t n = matrix.getNumCols();
      std::string name = "(" + std::to_string(s) + " * "  + matrix.getName() + ")";
      Matrix<T> l(name,m,n,0.0);
      for (uint32_t i = 0; i < m*n; i++) {
          l(i) = matrix(i) * s;
      }
      return l;
    }
    friend Matrix<T> operator/(T s, const Matrix<T>& matrix)
    {
      uint32_t m = matrix.getNumRows();
      uint32_t n = matrix.getNumCols();
      std::string name = "(" + std::to_string(s) + " / "  + matrix.getName() + ")";
      Matrix<T> l(name,m,n,0.0);
      std::vector<T> mat(n*m);
      for (uint32_t i = 0; i < m*n; i++)
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
      l.setArray(n, mat);
      return l;
    }
    //  Multiplying a vector
    Vector<T> operator*(const Vector<T>& v);
    //  Access operators
    T& operator()(const uint32_t& i, const uint32_t& j);
    const T& operator()(const uint32_t& i, const uint32_t& j) const;
    //  Flattened access
    T& operator()(const uint32_t& i);
    const T& operator()(const uint32_t& i) const;

    //  Various methods
    void print();
    std::string summary();
    Matrix<T> transpose() const;
    void transpose_inplace(bool inplace=true);
    T trace();


    //  Linear algebra tools
    Matrix<T> permutationMatrix(int& n, int* pivot);
    //  TODO:
    uint32_t getRank();
    bool isInvertible();
    void findSingularValues();
    Matrix<T> inverse();
    Matrix<T> pseudoInverse();
    //
    std::tuple<Matrix<T>,Matrix<T>,Matrix<T>> LU();
    Matrix<T> getL(const Matrix<T>& perm);
    Matrix<T> getU(const Matrix<T>& perm);
    std::tuple<Matrix<T>,Matrix<T>> QR();
    std::tuple<Matrix<T>,Matrix<T>,Matrix<T>> SVD();
    std::vector<T> getSingularValues();



  private:
    //  _m is the number of rows, _n is the number of columns
    uint32_t _m, _n;
    std::vector<T> _array;
    //  possible name for the matrix
    std::string _name;
    //  container for singular values
    std::vector<T> _singular_values;
    //  assign rank to -1 at initilization.
    uint32_t _rank = -1;
  };

  //  Various matrices
  template<typename T>
  Matrix<T> identity(uint32_t m);
  template<typename T>
  Matrix<T> zeroes(uint32_t m);
  template<typename T>
  Matrix<T> zeroes(uint32_t m, uint32_t n);
  template<typename T>
  Matrix<T> ones(uint32_t m);
  template<typename T>
  Matrix<T> ones(uint32_t m, uint32_t n);

  template<typename T>
  std::ostream& operator<<(std::ostream& os, const Matrix<T>& matrix);

  template class Matrix<double>;
  //template class Matrix<std::complex<double>>;

  //  Level 2 BLAS double
  //  DGEMV ()
  Vector<double> DGEMV(double& alpha, Matrix<double>& A,
      Vector<double>& x);
  Vector<double> DGEMV(double& alpha, Matrix<double>& A,
      Vector<double>& x, double& beta, Vector<double>& y);
}
