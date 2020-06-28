//------------------------------------------------------------------------------
//  matrix.h
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
  //--------------------------------------------------------------------------
  //  This matrix class acts as a general container for (m x n) matrices
  //  which are row ordered, i.e. m(i,j) picks the element of the ith row
  //  and the jth column.  It wraps several methods from BLAS and LAPACK.
  //--------------------------------------------------------------------------
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
    std::vector<T>* accessArray();
    //  get access to the beginning of the array
    T* data();
    std::vector<T> getRow(uint32_t i);
    std::vector<T> getCol(uint32_t i);
    std::vector<T> getSingularValues();
    std::string getInfo();
    int getFlag();
    uint32_t getRank();


    //  Setters
    void setName(std::string name);
    void setRow(uint32_t i, std::vector<T> row);
    void setCol(uint32_t i, std::vector<T> col);
    void setArray(uint32_t m, std::vector<T> mat);
    void setArray(std::vector<std::vector<T> > mat);
    void setSingularValues(std::vector<T> singular);
    void setInfo(std::string info);
    void setFlag(int flag);
    void setRank(uint32_t rank);

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
    //  TODO:
    bool isInvertible();
    void findSingularValues();
    Matrix<T> inverse();
    Matrix<T> pseudoInverse();
    std::tuple<Matrix<T>,Matrix<T>,Matrix<T>> LU();
    Matrix<T> getL(const Matrix<T>& perm);
    Matrix<T> getU(const Matrix<T>& perm);
    std::tuple<Matrix<T>,Matrix<T>> QR();
    std::tuple<Matrix<T>,Matrix<T>,Matrix<T>> SVD();

  private:
    //  _m is the number of rows, _n is the number of columns
    uint32_t _m, _n;
    std::vector<T> _array;
    //  possible name for the matrix
    std::string _name;
    //  container for singular values
    std::vector<T> _singular_values;
    //  conatiner for message status
    int _flag;
    //  container for messages
    std::string _info;
    //  assign rank to -1 at initilization.
    uint32_t _rank = -1;
  };

  //----------------------------------------------------------------------------
  //  Various Initializers
  //----------------------------------------------------------------------------
  Matrix<double> identity_d(uint32_t m);
  Matrix<double> zeros_d(uint32_t m);
  Matrix<double> zeros_d(uint32_t m, uint32_t n);
  Matrix<double> ones_d(uint32_t m);
  Matrix<double> ones_d(uint32_t m, uint32_t n);
  Matrix<double> permutationMatrix_d(const uint32_t& m,
                              const std::vector<uint32_t> pivot);
  template<typename T>
  std::ostream& operator<<(std::ostream& os, const Matrix<T>& matrix);

  template class Matrix<double>;
  //template class Matrix<std::complex<double>>;

  //----------------------------------------------------------------------------
  //  Level 2 BLAS methods
  //----------------------------------------------------------------------------
  Vector<double> DGEMV(double& alpha, Matrix<double>& A,
                       Vector<double>& x);
  void DGEMV(double& alpha, Matrix<double>& A, Vector<double>& x,
             double& beta, Vector<double>& y);
  Matrix<double> DGER(double& alpha, Vector<double>& x,
                      Vector<double>& y);
  void DGER(double& alpha, Vector<double>& x,
            Vector<double>& y, Matrix<double>& m);
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Level 3 BLAS methods
  //----------------------------------------------------------------------------
  Matrix<double> DGEMM(const double& alpha, const Matrix<double>& A,
                       const Matrix<double>& B);
  Matrix<double> DGEMM(const Matrix<double>& A,
                       const Matrix<double>& B);
  void DGEMM(const double& alpha, const Matrix<double>& A,
             const Matrix<double>& B, const double& beta,
             Matrix<double>& C);
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  LAPACK methods
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Linear least squares
  //----------------------------------------------------------------------------
  Matrix<double> DGELS(const Matrix<double>& A,
                       const Matrix<double>& B);
  Vector<double> DGELS(const Matrix<double>& A,
                       const Vector<double>& u);
  Matrix<double> DGELSY(const Matrix<double>& A,
                        const Matrix<double>& B);
  Vector<double> DGELSY(const Matrix<double>& A,
                        const Vector<double>& u);
  Matrix<double> DGELSD(const Matrix<double>& A,
                        const Matrix<double>& B);
  Vector<double> DGELSD(const Matrix<double>& A,
                        const Vector<double>& u);
  Matrix<double> DGELSS(const Matrix<double>& A,
                        const Matrix<double>& B);
  Vector<double> DGELSS(const Matrix<double>& A,
                        const Vector<double>& u);
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  LU decomposition
  //----------------------------------------------------------------------------
  std::vector<uint32_t> DGETRF(const Matrix<double>& A);
  Matrix<double> DGETRF_L_U(const Matrix<double>& A);
  std::tuple<Matrix<double>,Matrix<double>> DGETRF_LU(const Matrix<double>& A);
  std::tuple<Matrix<double>,Matrix<double>,Matrix<double>>
  DGETRF_PLU(const Matrix<double>& A);
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  QR decomposition
  //----------------------------------------------------------------------------
  Vector<double> DGEQRF(const Matrix<double>& A);
  Matrix<double> DORGQR(const Matrix<double>& A,
                        const Vector<double>& ref);
  std::tuple<Matrix<double>,Matrix<double>>
  DGEQRF_QR(const Matrix<double>& A);
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  SVD decomposition
  //----------------------------------------------------------------------------
  Vector<double> DGESVD(const Matrix<double>& A);
  std::tuple<Matrix<double>,Matrix<double>,Matrix<double>>
  DGESVD_SVD(const Matrix<double>& A);
  Vector<double> DGESDD(const Matrix<double>& A);
  std::tuple<Matrix<double>,Matrix<double>,Matrix<double>>
  DGESDD_SVD(const Matrix<double>& A);
  //----------------------------------------------------------------------------

}
