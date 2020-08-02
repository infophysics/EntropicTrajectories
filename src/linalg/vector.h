//------------------------------------------------------------------------------
// t_vector.h
// The Entropic Trajectories Framework
// -----------------------------------
// Copyright (C) [2020] by [N. Carrara]
// [ncarrara@albany.edu]
//
// Permission to use, copy, modify, and/or distribute this software for any
// purpose with or without fee is hereby granted.
//
// THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
// WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
// SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
// WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
// ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR
// IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
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
#include <iterator>

#include "utilities.h"
#include "log.h"

namespace ET
{
  /*! \class Vector
  *   The Vector class acts as a general container for (n)-dimensional
  *   Vectors.  It wraps several methods from BLAS level one such as
  *   ET::DSWAP(), ET::DSCAL(), ET::DCOPY, ET::DAXPY, ET::DNRM2, ET::DDOT,
  *   ET::DASUM, ET::IDAMAX() and ET::IDAMIN().
  */
  template<typename T>
  class Vector
  {
  public:
    //! Default constructor
    /*! Default constructor for a Vector object.
    */
    Vector();
    //! Destructor
    ~Vector();
    //! Constructor
    /*!
      @param t_vector a const Vector<T>& reference.
    */
    Vector(const Vector<T>& t_vector);
    //! Constructor
    /*!
      @param t_dim a size_t for the dimension of the vector space.
    */
    Vector(size_t t_dim);
    //! Constructor
    /*!
      @param t_name a std::string for the name of the Vector.
      @param t_dim a size_t for the dimension of the vector space.
    */
    Vector(std::string t_name, size_t t_dim);
    //! Constructor
    /*!
      @param t_vec an std::vector<T> containing the coefficients of the Vector.
    */
    Vector(std::vector<T> t_vec);
    //! Constructor
    /*!
      @param t_name a std::string for the name of the Vector.
      @param t_dim a size_t for the dimension of the vector space.
    */
    Vector(std::string t_name, std::vector<T> t_vec);
    //! Constructor
    /*!
      @param t_dim a size_t for the dimension of the vector space.
      @param t_init an initial value to set all of the elements to.
    */
    Vector(size_t t_dim, const T& t_init);
    //! Constructor
    /*!
      @param t_name a std::string for the name of thebrewathe sphinx Vector.
      @param t_dim a size_t for the dimension of the vector space.
      @param t_init an initial value to set all of the elements to.
    */
    Vector(std::string t_name, size_t t_dim, const T& t_init);
    //! Get dimension
    /*! get dimension.
        @return The dimension of the Vector.
    */
    size_t getDim() const;
    //! Get vec container
    /*! get vec container.
        @return A const reference of the std::vector<T> vec.
    */
    std::vector<T> getVec() const;
    //! Access vec container
    /*! access vec container.
        @return A raw pointer to the std::vector<T> container.
    */
    std::vector<T>* accessVec();
    //! Get pointer to the beginning of vecbrewathe sphinx
    /*! get vec pointer.  There are several ways we can interact with the
        m_vec attribute.  Some methods, such as Level one BLAS,
        will require access to the pointer for the first entry in m_vec.
        The data() method does just this.
        @returns A pointer to the beginning of the std::vector<T> vec.
    */
    T* data();
    //! Get name
    /*! get name.
        @return The name of the Vector.
    */
    std::string getName() const;
    //! Get flag
    /*! get flag.
     *  @return An int that classifies the type of information stored in
     *  info.
     */
    int getFlag() const;
    //! Get info
    /*! get info.
     *  @return The std::string info that contains relevant information.
     */
    std::string getInfo() const;
    //  Setters
    //! Set dimension
    /*! set dimension.  Sets the dimension \f$n\f$ of the the vector space
        \f$\mathbb{R}^n\f$
        @params t_dim a size_t for the dimension of the vector space.
    */
    void setDim(size_t t_dim);
    //! Set vec
    /*! set vec.  Sets the coefficients of the Vector.
        @param t_vec an std::vector<T> containing the coefficients of the Vector.
    */
    void setVec(std::vector<T> t_vec);
    //! Set name
    /*! set name.  Sets the name of the Vector.
        @param t_name a std::string for the name of the Vector.
    */
    void setName(std::string t_name);
    //! Set flag
    /*! set flag.  Sets the flag pertaining to info.
        @param t_flag an int that classifies the type of information stored
        in m_info.
    */
    void setFlag(int t_flag);
    //! Set info
    /*! set info.  Sets useful information pertaining to Vector.
        @param t_info an std::string containing useful messages.
    */
    void setInfo(std::string t_info);
    //  Operator overloads
    /*!
        @param t_vector A const Vector<T>& reference.
        @return A copy of the Vector
    */
    Vector<T>& operator=(const Vector<T>& t_vector);
    /*! Is equal.  Determines whether two vectors are equivalent.
        @param t_vector A const Vector<T>& reference.
        @return A boolean quantifying whether this Vector and vector
        are equivalent.
    */
    bool operator==(const Vector<T>& t_vector) const;
    /*! Is not equal.  Determines whether two vectors are not equal.
        @param t_vector A const Vector<T>& reference.
        @return A boolean quantifying whether this Vector and vector
        are not equivalent.
    */
    bool operator!=(const Vector<T>& t_vector) const;
    /*! Minus.  Constructs a copy of the vector \f$\vec{v}\f$ multiplied by
        minus one, \f$-\vec{v}\f$.
        @return A copy of this Vector multiplied by \f$-1\f$.
    */
    Vector<T> operator-() const;
    /*! Sum.  Adds two vectors together.
        @param t_vector A const Vector<T>& reference.
        @return The sum of this Vector and vector.
    */
    Vector<T> operator+(const Vector<T>& t_vector) const;
    /*! Sum equals.  Adds a vector to the current one.
        @param t_vector A const Vector<T>& reference.
        @return This Vector plus vector
    */
    Vector<T>& operator+=(const Vector<T>& t_vector);
    /*! Difference.  Subtracts two vectors.
        @param t_vector A const Vector<T>& reference.
        @return The difference of this Vector and vector.
    */
    Vector<T> operator-(const Vector<T>& t_vector) const;
    /*! Difference equals.  Subtracts a vector from the current one.
        @param t_vector A const Vector<T>& reference.
        @return This vector minus vector.
    */
    Vector<T>& operator-=(const Vector<T>& t_vector);
    //  Scalar product
    /*! Scalar product.  Computes the scalar product between two vectors,
        \f$\beta = \vec{v}\cdot \vec{u}\f$.
        @param t_vector A const Vector<T>& reference.
        @return The scalar product between this Vector and vector.
    */
    T operator*(const Vector<T>& t_vector) const;
    /*! Scalar product.  Computes the scalar product between two vectors,
        \f$\beta = \vec{v}\cdot\vec{u}\f$.
        @param t_vector A const Vector<T>& reference.
        @return The scalar product between this Vector and vector.
    */
    T dot(const Vector<T>& t_vector) const;
    //  Scalar operators
    /*! Sum of a scalar.  Adds a scalar times an identity vector
        to another vector, \f$\vec{v} = \vec{u} +
        \alpha \mathbf{1}\f$.
        @param t_s A const T& reference.
        @return This Vector plus the scalar.
    */
    Vector<T> operator+(const T& t_s) const;
    /*! Difference of a scalar.  Subtracts a scalar times an identity vector
        to another vector, \f$\vec{v} = \vec{u} -
        \alpha \mathbf{1}\f$.
        @param t_s A const T& reference.
        @return This Vector minus the scalar.
    */
    Vector<T> operator-(const T& t_s) const;
    /*! Scalar multiplication.  Multiplies a scalar with a vector,
        \f$\vec{v} = \alpha\vec{u}\f$.
        @param t_s A const T& reference.
        @return This Vector times the scalar.
    */
    Vector<T> operator*(const T& t_s) const;
    /*! Scalar division.  Divides a vector by a scalar,
        \f$\vec{v} = \vec{u}/\alpha
        = (\alpha)^{-1}\vec{u}\f$.
        @param t_s A const T& reference.
        @return This Vector divided by the scalar.
    */
    Vector<T> operator/(const T& t_s) const;
    /*! Sum equals of a scalar.  Adds a scalar times an identity vector
        to this vector, \f$\vec{v} \rightarrow \vec{v}' = \vec{v} +
        \alpha \mathbf{1}\f$.
        @param t_s A const T& reference.
        @return This Vector plus the scalar.
    */
    Vector<T>& operator+=(const T& t_s);
    /*! Difference equals of a scalar.  Subtracts a scalar times an identity vector
        to this vector, \f$\vec{v} \rightarrow \vec{v}' = \vec{v} -
        \alpha \mathbf{1}\f$.
        @param t_s A const T& reference.
        @return This Vector minus the scalar.
    */
    Vector<T>& operator-=(const T& t_s);
    /*! Scalar equals multiplication.  Multiplies a scalar with this vector,
        \f$\vec{v} \rightarrow \vec{v}' = \alpha\vec{v}\f$.
        @param t_s A const T& reference.
        @return This Vector times the scalar.
    */
    Vector<T>& operator*=(const T& t_s);
    /*! Scalar equals division.  Divides this vector by a scalar,
        \f$\vec{v} \rightarrow \vec{v}' = \vec{v}/\alpha
        = (\alpha)^{-1}\vec{v}\f$.
        @param t_s A const T& reference.
        @return This Vector divided by the scalar.
    */
    Vector<T>& operator/=(const T& t_s);
    //  Overloads of scalar operations from the left.
    //  since we are trying to friend a template argument,
    //  the friend method must be defined within the class block.
    /*! Sum of a scalar from the left.  Adds a scalar times an identity vector
        to another vector, \f$\vec{v} = \alpha \mathbf{1} + \vec{u}\f$.
        @param t_s A const T& reference.
        @param t_vector A const Vector<T>& reference.
        @return This Vector plus the scalar.
    */
    friend Vector<T> operator+(T t_s, const Vector<T>& t_vector)
    {
      size_t dim = t_vector.getDim();
      std::string name = "(" + std::to_string(t_s) + " + "
                         + t_vector.getName() + ")";
      Vector<T> v(name,dim,0.0);
      for (auto i = 0; i < dim; i++) {
          v(i) = t_vector(i) + t_s;
      }
      return v;
    }
    /*! Difference of a scalar from the left.
        Subtracts a vector from a scalar times an identity vector,
        \f$\vec{v} = \alpha \mathbf{1} - \vec{u}\f$.
        @param t_s A const T& reference.
        @param t_vector A const Vector<T>& reference.
        @return A scalar minus vector.
    */
    friend Vector<T> operator-(T t_s, const Vector<T>& t_vector)
    {
      size_t dim = t_vector.getDim();
      std::string name = "(" + std::to_string(t_s) + " - "
                       + t_vector.getName() + ")";
      Vector<T> v(name,dim,0.0);
      for (auto i = 0; i < dim; i++) {
          v(i) = t_s - t_vector(i);
      }
      return v;
    }
    /*! Scalar multiplication from the left.  Multiplies a scalar with a vector,
        \f$\vec{v} = \alpha\vec{u}\f$.
        @param t_s A const T& reference.
        @return This Vector times the scalar.
    */
    friend Vector<T> operator*(T t_s, const Vector<T>& t_vector)
    {
      size_t dim = t_vector.getDim();
      std::string name = "(" + std::to_string(t_s) + " * "
                       + t_vector.getName() + ")";
      Vector<T> v(name,dim,0.0);
      for (auto i = 0; i < dim; i++) {
          v(i) = t_vector(i) * t_s;
      }
      return v;
    }
    /*! Scalar division from the left.  Divides each element of the vector
        by a scalar,
        \f$\forall v_i \in \vec{v} : v_i \rightarrow v_i' = \alpha/v_i\f$.
        @param t_s A const T& reference.
        @return A vector whose components are the scalar divided by
        the components of this Vector.
    */
    friend Vector<T> operator/(T t_s, const Vector<T>& t_vector)
    {
      size_t dim = t_vector.getDim();
      std::string name = "(" + std::to_string(t_s) + " / "
                       + t_vector.getName() + ")";
      Vector<T> v(name,dim,0.0);
      std::vector<T> vec(dim);
      for (auto i = 0; i < dim; i++) {
        if (t_vector(i) == 0) {
          return v;
        }
        else {
         vec[i] = t_s / t_vector(i);
        }
      }
      v.setVec(vec);
      return v;
    }

    //  Access operators
    /*! Access operator.  An access operator for changing entries of the vector.
        @param t_i A const size_t& reference for the index of
        the desired coefficient.
        @return The ith-coefficient.
    */
    T& operator()(const size_t& t_i);
    /*! Access operator.  An access operator for retrieving coefficients
        of the vector.
        @param t_i A const size_t& reference for the index of
        the desired coefficient.
        @return The ith-coefficient.
    */
    const T& operator()(const size_t& t_i) const;
    /*! Summary.  Returns a summary of information about the Vector
        @return An std::string of information about the Vector.
     */
    const std::string summary();
    //! Remove value
    /*! Remove a value to decrease the dimension of the vector.
     *  @param t_index The location to remove the value.
     */
    Vector<T> removeVal(const size_t t_index) const;
    //! Remove value inplace
    /*! Remove a value to decrease the dimension of the vector.
     *  @param t_index The location to remove the value.
     */
    void removeValInplace(const size_t t_index);
    //! Add Value
    /*! Add a new value to increase the dimension of the vector.
     *  @param t_index The location to insert the value.
     *  @param t_value The value to insert.
     */
    Vector<T> addVal(const size_t t_index, const T t_value) const;
    //! Add Value inplace
    /*! Add a new value to increase the dimension of the vector.
     *  @param t_index The location to insert the value.
     *  @param t_value The value to insert.
     */
    void addValInplace(const size_t t_index, const T t_value);

  private:
    /*!  Dimension \f$n\f$ of the vector space \f$\mathbb{R}^n\f$.
         Defaulted to 0.
     */
    size_t m_dim {0};
    /*!  Container for the coefficients in \f$R^n\f$.
         Defaulted to a single entry whose value is 0.
     */
    std::vector<T> m_vec {{0}};
    /*!  Name of the Vector.
         Defaulted to an empty string.
     */
    std::string m_name {""};
    /*!  A flag for the type of message stored in m_info.
         Defaulted to 0.
    */
    int m_flag {0};
    /*!  Useful runtime information.
         Defaulted to an empty string.
     */
    std::string m_info {""};

  public:
    //  Inheriting iterators from std::vector
    using iterator = typename std::vector<T>::iterator;
    using const_iterator = typename std::vector<T>::const_iterator;

    iterator begin()              { return m_vec.begin(); }
    iterator end()                { return m_vec.end(); }
    const_iterator cbegin() const { return m_vec.cbegin(); }
    const_iterator cend() const   { return m_vec.cend(); }
  };

  //----------------------------------------------------------------------------
  //  Special Vector initializers
  //----------------------------------------------------------------------------
  template<typename T>
  Vector<T> zeroes(size_t t_dim);
  template<typename T>
  Vector<T> ones(size_t t_dim);
  //----------------------------------------------------------------------------

  template class Vector<double>;

  //----------------------------------------------------------------------------
  //  Level 1 BLAS methods
  //----------------------------------------------------------------------------
  /*! DSWAP. Swap the contents of two Vectors using the BLAS level one
   *  routine **cblas_dswap**
   *  from the Intel mkl library.
   *  @param v A Vector<double>& reference
   *  @param u A Vector<double>& reference
   *  @return void
   */
  void DSWAP(Vector<double>& t_v, Vector<double>& t_u);
  /*! DSCAL.  Multiplies a vector by a scalar using the BLAS level one
   *  routine **cblas_dscal**
   *  from the Intel mkl library.
   *  @param v A Vector<double>& reference
   *  @param scale A const double& reference
   *  @return void
   */
  void DSCAL(Vector<double>& t_v, const double& t_scale);
  /*! DCOPY.  Copies the contents of a vector into a new vector
   *  using the BLAS level one routine **cblas_dcopy**
   *  from the Intel mkl library.
   *  @param v A Vector<double>& reference
   *  @return A Copy of the vector v.
   */
  Vector<double> DCOPY(Vector<double>& t_v);
  /*! DCOPY.  Copies the contents of a vector into another
   *  using the BLAS level one routine **cblas_dcopy**
   *  from the Intel mkl library.
   *  @param v A Vector<double>& reference which is the vector to be copied.
   *  @param u A Vector<double>& reference which is the vector to be replaced
   *  with the contents of v.
   *  @return void
   */
  void DCOPY(Vector<double>& t_v, Vector<double>& t_u);
  /*! DAXPY.  Copies the contents of a vector times a scalar into another
   *  using the BLAS level one routine **cblas_daxpy**
   *  from the Intel mkl library.
   *  @param v A Vector<double>& reference which is the vector to be copied.
   *  @param scale A const double& reference
   *  @param u A Vector<double>& reference which is the vector to be replaced
   *  with the contents of v times the scalar.
   *  @return void
   */
  void DAXPY(Vector<double>& t_v, const double& scale, Vector<double>& t_u);
  /*! DDOT.  Computes the dot product between two vectors
   *  using the BLAS level one routine **cblas_ddot**
   *  from the Intel mkl library.
   *  @param v A Vector<double>& reference
   *  @param u A Vector<double>& reference
   *  @return The dot product of v and u.
   */
  double DDOT(Vector<double>& t_v, Vector<double>& t_u);
  /*! DNRM2.  Computes the Euclidean norm of a Vector
   *  using the BLAS level one routine **cblas_dnrm2**
   *  from the Intel mkl library.
   *  @param v A Vector<double>& reference
   *  @return The Euclidean norm of v.
   */
  double DNRM2(Vector<double>& t_v);
  /*! DASUM.  Computes the sum of the absolute value of each of the coefficients
   *  of a vector, \f$\alpha = \sum_{i=1}^n |v_i|\f$,
   *  using the BLAS level one routine **cblas_dasum**
   *  from the Intel mkl library.
   *  @param v A Vector<double>& reference
   *  @return The sum of the absolute value of the elements of v.
   */
  double DASUM(Vector<double>& t_v);
  /*! IDAMAX.  Returns the index of the coefficient in a vector with the
   *  largest value using the BLAS level one routine **cblas_idamax**
   *  from the Intel mkl library.
   *  @param v A Vector<double>& reference
   *  @return The index of the largest value of v.
   */
  size_t IDAMAX(Vector<double>& t_v);
  /*! IDAMIN.  Returns the index of the coefficient in a vector with the
   *  smallest value using the BLAS level one routine **cblas_idamax**
   *  from the Intel mkl library.
   *  @param v A Vector<double>& reference
   *  @return The index of the smallest value of v.
   */
  size_t IDAMIN(Vector<double>& t_v);
  //----------------------------------------------------------------------------
}
