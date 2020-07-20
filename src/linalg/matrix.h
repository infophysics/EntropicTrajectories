//------------------------------------------------------------------------------
//  matrix.h
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
  //! Matrix Class
  /*! This matrix class acts as a wrapper for several BLAS and LAPACK
   *  routines that will be necessary in the ET framework.
   */
  template<typename T>
  class Matrix
  {
  public:
    //! Default Constructor
    /*! Default constructor for a matrix object.
     */
    Matrix();
    //! Destructor
    /*  Destructor for a matrix object.
     */
    ~Matrix();
    //! Copy constructor
    /*  Override for the copy constructor.
     */
    Matrix(const Matrix<T>& t_matrix);
    //! Constructor
    /*!
     *
     */
    Matrix(size_t t_m);
    //! Constructor
    /*!
     *
     */
    Matrix(std::string t_name, size_t t_m);
    //! Constructor
    /*!
     *
     */
    Matrix(size_t t_m, size_t t_n);
    //! Constructor
    /*!
     *
     */
    Matrix(std::string t_name, size_t t_m, size_t t_n);
    //! Constructor
    /*!
     *
     */
    Matrix(size_t t_m, size_t t_n, const T& t_init);
    //! Constructor
    /*!
     *
     */
    Matrix(std::string t_name, size_t t_m, size_t t_n, const T& t_init);
    //! Constructor
    /*!
     *
     */
    Matrix(size_t t_m, std::vector<T> t_flat);
    //! Constructor
    /*!
     *
     */
    Matrix(std::string t_name, size_t t_m, std::vector<T> t_flat);
    //! Constructor
    /*!
     *
     */
    Matrix(size_t t_m, size_t t_n, std::vector<T> t_flat);
    //! Constructor
    /*!
     *
     */
    Matrix(std::string t_name, size_t t_m,
      size_t t_n, std::vector<T> t_flat);
      //! Constructor
      /*!
       *
       */
    Matrix(std::string t_name, size_t t_m, size_t t_n, T* t_array);
    //! Constructor
    /*!
     *
     */
    Matrix(std::vector<std::vector<T>> t_array);
    //! Constructor
    /*!
     *
     */
    Matrix(std::string t_name, std::vector<std::vector<T>> t_array);

    //  Getters
    size_t getNumRows() const;
    size_t getNumCols() const;
    std::string getName() const;
    //  get a const reference to array
    std::vector<T> getArray() const;
    //  get access to array
    std::vector<T>* accessArray();
    //  get access to the beginning of the array
    T* data();
    std::vector<T> getRow(size_t t_i);
    std::vector<T> getCol(size_t t_i);
    std::vector<T> getSingularValues();
    std::string getInfo();
    int getFlag();
    size_t getRank();


    //  Setters
    void setName(std::string t_name);
    void setRow(size_t t_i, std::vector<T> t_row);
    void setCol(size_t t_i, std::vector<T> t_col);
    void setArray(size_t t_m, std::vector<T> t_mat);
    void setArray(std::vector<std::vector<T>> t_mat);
    void setSingularValues(std::vector<T> t_singular);
    void setInfo(std::string t_info);
    void setFlag(int t_flag);
    void setRank(size_t t_rank);

    //  Operator overloads
    Matrix<T>& operator=(const Matrix<T>& t_matrix);
    bool operator==(const Matrix<T>& t_matrix) const;
    bool operator!=(const Matrix<T>& t_matrix) const;
    Matrix<T> operator-() const;
    //  Matrix algebra
    Matrix<T> operator+(const Matrix<T>& t_matrix) const;
    Matrix<T>& operator+=(const Matrix<T>& t_matrix);
    Matrix<T> operator-(const Matrix<T>& t_matrix) const;
    Matrix<T>& operator-=(const Matrix<T>& t_matrix);
    //  Using CBLAS for multiplication
    Matrix<T> operator*(const Matrix<T>& t_matrix) const;
    Matrix<T> brutem_mul(const Matrix<T>& t_matrix) const;
    Matrix<T>& operator*=(const Matrix<T>& t_matrix);
    //  Scalar algebra
    Matrix<T> operator+(const T& t_s) const;
    Matrix<T> operator-(const T& t_s) const;
    Matrix<T> operator*(const T& t_s) const;
    Matrix<T> operator/(const T& t_s) const;
    Matrix<T>& operator+=(const T& t_s);
    Matrix<T>& operator-=(const T& t_s);
    Matrix<T>& operator*=(const T& t_s);
    Matrix<T>& operator/=(const T& t_s);
    //  Overloads of scalar operations from the left.
    //  since we are trying to friend a template argument,
    //  the friend method must be defined within the class block.
    friend Matrix<T> operator+(T t_s, const Matrix<T>& t_matrix)
    {
      size_t m = t_matrix.getNumRows();
      size_t n = t_matrix.getNumCols();
      std::string name = "(" + std::to_string(t_s) + "I + "
                         + t_matrix.getName() + ")";
      Matrix<T> l(name,m,n,0.0);
      for (size_t i = 0; i < m*n; i++) {
          l(i) = t_matrix(i) + t_s;
      }
      return l;
    }
    friend Matrix<T> operator-(T t_s, const Matrix<T>& t_matrix)
    {
      size_t m = t_matrix.getNumRows();
      size_t n = t_matrix.getNumCols();
      std::string name = "(" + std::to_string(t_s) + "I - "
                         + t_matrix.getName() + ")";
      Matrix<T> l(name,m,n,0.0);
      for (size_t i = 0; i < m*n; i++) {
          l(i) = t_s - t_matrix(i);
      }
      return l;
    }
    friend Matrix<T> operator*(T t_s, const Matrix<T>& t_matrix)
    {
      size_t m = t_matrix.getNumRows();
      size_t n = t_matrix.getNumCols();
      std::string name = "(" + std::to_string(t_s) + " * "
                         + t_matrix.getName() + ")";
      Matrix<T> l(name,m,n,0.0);
      for (size_t i = 0; i < m*n; i++) {
          l(i) = t_matrix(i) * t_s;
      }
      return l;
    }
    friend Matrix<T> operator/(T t_s, const Matrix<T>& t_matrix)
    {
      size_t m = t_matrix.getNumRows();
      size_t n = t_matrix.getNumCols();
      std::string name = "(" + std::to_string(t_s) + " / "
                         + t_matrix.getName() + ")";
      Matrix<T> l(name,m,n,0.0);
      std::vector<T> mat(n*m);
      for (size_t i = 0; i < m*n; i++)
      {
        if (t_matrix(i) == 0) {
          return l;
        }
        else {
          mat[i] = t_s / t_matrix(i);
        }
      }
      l.setArray(n, mat);
      return l;
    }
    //  Multiplying a vector
    Vector<T> operator*(const Vector<T>& t_v);
    //  Access operators
    T& operator()(const size_t& t_i, const size_t& t_j);
    const T& operator()(const size_t& t_i, const size_t& t_j) const;
    //  Flattened access
    T& operator()(const size_t& t_i);
    const T& operator()(const size_t& t_i) const;

    //  Various methods
    void print();
    const std::string summary();
    Matrix<T> transpose() const;
    void transpose_inplace(bool t_inplace=true);
    T trace();

    //  Linear algebra tools
    //  TODO:
    // bool isInvertible();
    // void findSingularValues();
    // Matrix<T> inverse();
    // Matrix<T> pseudoInverse();
    // std::tuple<Matrix<T>,Matrix<T>,Matrix<T>> LU();
    // Matrix<T> getL(const Matrix<T>& perm);
    // Matrix<T> getU(const Matrix<T>& perm);
    // std::tuple<Matrix<T>,Matrix<T>> QR();
    // std::tuple<Matrix<T>,Matrix<T>,Matrix<T>> SVD();

  private:
    //  m_m is the number of rows, m_n is the number of columns
    size_t m_m, m_n;
    std::vector<T> m_array;
    //  possible name for the t_matrix
    std::string m_name;
    //  container for singular values
    std::vector<T> m_singular_values;
    //  conatiner for message status
    int m_flag;
    //  container for messages
    std::string m_info;
    //  assign rank to -1 at initilization.
    size_t m_rank = -1;
  };

  //----------------------------------------------------------------------------
  //  Various Initializers
  //----------------------------------------------------------------------------
  Matrix<double> identity_d(size_t t_m);
  Matrix<double> zeros_d(size_t t_m);
  Matrix<double> zeros_d(size_t t_m, size_t t_n);
  Matrix<double> ones_d(size_t t_m);
  Matrix<double> ones_d(size_t t_m, size_t t_n);
  Matrix<double> permutationMatrix_d(const size_t& t_m,
                              const std::vector<size_t> t_pivot);
  template<typename T>
  std::ostream& operator<<(std::ostream& t_os, const Matrix<T>& t_matrix);

  template class Matrix<double>;
  //template class Matrix<std::complex<double>>;

  //----------------------------------------------------------------------------
  //  Level 2 BLAS methods
  //----------------------------------------------------------------------------
  Vector<double> DGEMV(double& t_alpha, Matrix<double>& t_A,
                       Vector<double>& t_x);
  void DGEMV(double& t_alpha, Matrix<double>& t_A, Vector<double>& t_x,
             double& t_beta, Vector<double>& t_y);
  Matrix<double> DGER(double& t_alpha, Vector<double>& t_x,
                      Vector<double>& t_y);
  void DGER(double& t_alpha, Vector<double>& t_x,
            Vector<double>& t_y, Matrix<double>& t_m);
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Level 3 BLAS methods
  //----------------------------------------------------------------------------
  Matrix<double> DGEMM(const double& t_alpha, const Matrix<double>& t_A,
                       const Matrix<double>& t_B);
  Matrix<double> DGEMM(const Matrix<double>& t_A,
                       const Matrix<double>& t_B);
  void DGEMM(const double& t_alpha, const Matrix<double>& t_A,
             const Matrix<double>& t_B, const double& t_beta,
             Matrix<double>& t_C);
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  LAPACK methods
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Linear least squares
  //----------------------------------------------------------------------------
  Matrix<double> DGELS(const Matrix<double>& t_A,
                       const Matrix<double>& t_B);
  Vector<double> DGELS(const Matrix<double>& t_A,
                       const Vector<double>& t_u);
  Matrix<double> DGELSY(const Matrix<double>& t_A,
                        const Matrix<double>& t_B);
  Vector<double> DGELSY(const Matrix<double>& t_A,
                        const Vector<double>& t_u);
  Matrix<double> DGELSD(const Matrix<double>& t_A,
                        const Matrix<double>& t_B);
  Vector<double> DGELSD(const Matrix<double>& t_A,
                        const Vector<double>& t_u);
  Matrix<double> DGELSS(const Matrix<double>& t_A,
                        const Matrix<double>& t_B);
  Vector<double> DGELSS(const Matrix<double>& t_A,
                        const Vector<double>& t_u);
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  LU decomposition
  //----------------------------------------------------------------------------
  std::vector<size_t> DGETRF(const Matrix<double>& t_A);
  Matrix<double> DGETRF_L_U(const Matrix<double>& t_A);
  std::tuple<Matrix<double>,Matrix<double>>
  DGETRF_LU(const Matrix<double>& t_A);
  std::tuple<Matrix<double>,Matrix<double>,Matrix<double>>
  DGETRF_PLU(const Matrix<double>& t_A);
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  QR decomposition
  //----------------------------------------------------------------------------
  Vector<double> DGEQRF(const Matrix<double>& t_A);
  Matrix<double> DORGQR(const Matrix<double>& t_A,
                        const Vector<double>& t_ref);
  std::tuple<Matrix<double>,Matrix<double>>
  DGEQRF_QR(const Matrix<double>& t_A);
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  SVD decomposition
  //----------------------------------------------------------------------------
  Vector<double> DGESVD(const Matrix<double>& t_A);
  std::tuple<Matrix<double>,Matrix<double>,Matrix<double>>
  DGESVD_SVD(const Matrix<double>& t_A);
  Vector<double> DGESDD(const Matrix<double>& t_A);
  std::tuple<Matrix<double>,Matrix<double>,Matrix<double>>
  DGESDD_SVD(const Matrix<double>& t_A);
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Inverses
  //----------------------------------------------------------------------------
  Matrix<double> DGETRI(const Matrix<double>& t_A);
  //----------------------------------------------------------------------------
}
