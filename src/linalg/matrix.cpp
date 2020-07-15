//------------------------------------------------------------------------------
//  matrix.cpp
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
#include "vector.h"
#include "matrix.h"

namespace ET
{
  //----------------------------------------------------------------------------
  //  Matrix constructors
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Default constructor
  //    sets name = " ", and m_m, m_n = 0
  //----------------------------------------------------------------------------
  template<typename T>
  Matrix<T>::Matrix() : m_name(" "), m_m(0), m_n(0)
  {
    //std::cout << "\nMatrix created at location " << this;
  }
  //----------------------------------------------------------------------------
  //  Default destructor
  //----------------------------------------------------------------------------
  template<typename T>
  Matrix<T>::~Matrix()
  {
    //std::cout << "\nMatrix " + m_name + " at location " << this << " destroyed.";
  }
  //----------------------------------------------------------------------------
  //  Copy constructor
  //    Does not delete the copied object.
  //----------------------------------------------------------------------------
  template<typename T>
  Matrix<T>::Matrix(const Matrix<T>& t_matrix)
  {
    //std::cout << "\nMatrix created at location " << this;
    m_array = t_matrix.getArray();
    m_m = t_matrix.getNumRows();
    m_n = t_matrix.getNumCols();
    m_name = t_matrix.getName();
  }
  //----------------------------------------------------------------------------
  //  Constructors with various sets of arguments, such as,
  //    std::string                 name,
  //    size_t                    m,
  //    size_t                    n,
  //    std::vector<T>              t_flatted t_array,
  //    std::vector<std::vector<T>> 2d t_array,
  //    T*                          C style t_flattened t_array,
  //    const T&                    t_initial value for all elements.
  //  If any sizes are specified, the internal t_array element '_array'
  //  is resized accordingly, otherwise it is left uninstantiated.
  //----------------------------------------------------------------------------
  template<typename T>
  Matrix<T>::Matrix(size_t t_m) : m_m(t_m), m_n(t_m), m_name(" ")
  {
    //std::cout << "\nMatrix created at location " << this;
    m_array.resize(m_m*m_n,0.0);
  }
  template<typename T>
  Matrix<T>::Matrix(std::string t_name, size_t t_m)
  : m_m(t_m), m_n(t_m), m_name(t_name)
  {
    //std::cout << "\nMatrix created at location " << this;
    m_array.resize(m_m*m_n,0.0);
  }
  template<typename T>
  Matrix<T>::Matrix(size_t t_m, size_t t_n) : m_m(t_m), m_n(t_n), m_name(" ")
  {
    //std::cout << "\nMatrix created at location " << this;
    m_array.resize(m_m*m_n,0.0);
  }
  template<typename T>
  Matrix<T>::Matrix(std::string t_name, size_t t_m, size_t t_n)
  : m_m(t_m), m_n(t_n), m_name(t_name)
  {
    //std::cout << "\nMatrix created at location " << this;
    m_array.resize(m_m*m_n,0.0);
  }
  template<typename T>
  Matrix<T>::Matrix(size_t t_m, size_t t_n, const T& t_init)
  : m_m(t_m), m_n(t_n), m_name(" ")
  {
    //std::cout << "\nMatrix created at location " << this;
    m_array.resize(m_m*m_n, t_init);
  }
  template<typename T>
  Matrix<T>::Matrix(std::string t_name, size_t t_m,
    size_t t_n, const T& t_init)
  : m_m(t_m), m_n(t_n), m_name(t_name)
  {
    //std::cout << "\nMatrix created at location " << this;
    m_array.resize(m_m*m_n, t_init);
  }
  //  Notice that the following methods do not MOVE the vectors so that they
  //  change ownership.  Instead they are copied into _array.
  template<typename T>
  Matrix<T>::Matrix(size_t t_m, std::vector<T> t_flat)
  : m_m(t_m), m_n(t_m), m_name(" "), m_array(t_flat)
  {
    //std::cout << "\nMatrix created at locatiot_n " << this;
  }
  template<typename T>
  Matrix<T>::Matrix(std::string t_name, size_t t_m, std::vector<T> t_flat)
  : m_m(t_m), m_n(t_m), m_name(t_name), m_array(t_flat)
  {
    //std::cout << "\nMatrix created at location " << this;
  }
  template<typename T>
  Matrix<T>::Matrix(size_t t_m, size_t t_n, std::vector<T> t_flat)
  : m_m(t_m), m_n(t_n), m_name(" "), m_array(t_flat)
  {
    //std::cout << "\nMatrix created at location " << this;
  }
  template<typename T>
  Matrix<T>::Matrix(std::string t_name, size_t t_m, size_t t_n, std::vector<T> t_flat)
  : m_m(t_m), m_n(t_n), m_name(t_name), m_array(t_flat)
  {
    //std::cout << "\nMatrix created at location " << this;
  }
  template<typename T>
  Matrix<T>::Matrix(std::string t_name, size_t t_m, size_t t_n, T* t_array)
  : m_m(t_m), m_n(t_n), m_name(t_name)
  {
    //std::cout << "\nMatrix created at location " << this;
    std::vector<T> flat(t_array, t_array + m_m*m_n);
    m_array = flat;
  }
  template<typename T>
  Matrix<T>::Matrix(std::vector<std::vector<T>> t_array)
  : m_m(t_array.size()), m_n(t_array[0].size()), m_name(" ")
  {
    //std::cout << "\nMatrix created at location " << this;
    std::vector<T> flat;
    for (auto i = 0; i < m_m; i++) {
      flat.insert(end(flat),begin(t_array[i]),end(t_array[i]));
    }
    m_array = flat;
  }
  template<typename T>
  Matrix<T>::Matrix(std::string t_name, std::vector<std::vector<T>> t_array)
  : m_m(t_array.size()), m_n(t_array[0].size()), m_name(t_name)
  {
    //std::cout << "\nMatrix created at location " << this;
    std::vector<std::pair<size_t,size_t>> rows;
    bool well_defined = checkConsistency(t_array,rows);
    if (well_defined == false) {
      setInfo(MATRIX_INCONSISTENT_ARRAY(rows));
    }
    std::vector<T> flat;
    for (auto i = 0; i < m_m; i++) {
      flat.insert(end(flat),begin(t_array[i]),end(t_array[i]));
    }
    m_array = flat;
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Getters and Setters
  //  Each attribute comes with its own setters and getters.  There are a few
  //  additional ones such as 'getRow' and 'getCol' which each return an
  //  std::vector<T>.
  //----------------------------------------------------------------------------
  template<typename T>
  size_t Matrix<T>::getNumRows() const
  {
    return m_m;
  }
  template<typename T>
  size_t Matrix<T>::getNumCols() const
  {
    return m_n;
  }
  template<typename T>
  std::string Matrix<T>::getName() const
  {
    return m_name;
  }
  //  When 'getArray()' is called it will usually create a copy of _array.
  template<typename T>
  std::vector<T> Matrix<T>::getArray() const
  {
    return m_array;
  }
  //  In order to return the '_array' attribute so that it can be manipulated
  //  by other methods, such as those which utilize BLAS and LAPACK functions,
  //  we use 'accessArray()' to return a pointer to '_array'.
  template<typename T>
  std::vector<T>* Matrix<T>::accessArray()
  {
    return &m_array;
  }
  //  get access to the beginning of the t_array
  template<typename T>
  T* Matrix<T>::data()
  {
    return m_array.data();
  }
  template<typename T>
  std::vector<T> Matrix<T>::getRow(size_t t_i)
  {
    std::vector<T> row(m_m,0.0);
    //  check that row exists, if not log error and send zero vector
    if (t_i >= m_m) {
      m_info = MATRIX_OUT_OF_BOUNDS(0,m_n,t_i,m_name);
      return row;
    }
    for (auto j = 0; j < m_m; j++) {
      row[j] = m_array[t_i*m_m + j];
    }
    return row;
  }
  template<typename T>
  std::vector<T> Matrix<T>::getCol(size_t t_i)
  {
    std::vector<T> col(m_m,0.0);
    for (auto j = 0; j < m_n; j++) {
      col[j] = m_array[j*m_n + t_i];
    }
    return col;
  }
  template<typename T>
  std::vector<T> Matrix<T>::getSingularValues()
  {
    return m_singular_values;
  }
  template<typename T>
  std::string Matrix<T>::getInfo()
  {
    return m_info;
  }
  template<typename T>
  int Matrix<T>::getFlag()
  {
    return m_flag;
  }
  template<typename T>
  size_t Matrix<T>::getRank()
  {
    return m_rank;
  }
  //  Setters
  template<typename T>
  void Matrix<T>::setName(std::string t_name)
  {
    m_name = t_name;
  }
  template<typename T>
  void Matrix<T>::setRow(size_t t_i, std::vector<T> t_row)
  {
    for (auto j = 0; j < m_n; j++) {
      m_array[t_i*m_n + j] = t_row[j];
    }
  }
  template<typename T>
  void Matrix<T>::setCol(size_t t_i, std::vector<T> t_col)
  {
    for (auto j = 0; j < m_m; j++) {
      m_array[j*m_n + t_i] = t_col[j];
    }
  }
  template<typename T>
  void Matrix<T>::setArray(size_t t_m, std::vector<T> t_mat)
  {
    m_n = t_mat.size()/t_m;
    m_m = t_m;
    m_array = t_mat;
  }
  template<typename T>
  void Matrix<T>::setArray(std::vector<std::vector<T>> t_mat)
  {
    std::vector<std::pair<size_t,size_t>> rows;
    bool well_defined = checkConsistency(t_mat,rows);
    if (well_defined == false)
    {
      setInfo(MATRIX_INCONSISTENT_ARRAY(rows));
    }
    m_m = t_mat.size();
    m_n = t_mat[0].size();
    m_array.resize(m_m*m_n);
    for (auto i = 0; i < m_m; i++) {
      for (auto j = 0; j < m_n; j++) {
        m_array[i*m_n + j] = t_mat[i][j];
      }
    }
  }
  template<typename T>
  void Matrix<T>::setSingularValues(std::vector<T> t_singular)
  {
    m_singular_values = t_singular;
  }
  template<typename T>
  void Matrix<T>::setInfo(std::string t_info)
  {
    m_info = t_info;
  }
  template<typename T>
  void Matrix<T>::setFlag(int t_flag)
  {
    m_flag = t_flag;
  }
  template<typename T>
  void Matrix<T>::setRank(size_t t_rank)
  {
    m_rank = t_rank;
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Operator overloads
  //----------------------------------------------------------------------------
  template<typename T>
  Matrix<T>& Matrix<T>::operator=(const Matrix<T>& t_matrix)
  {
    if (&t_matrix == this) {
      return *this;
    }
    m_m = t_matrix.getNumRows();
    m_n = t_matrix.getNumCols();
    m_name = t_matrix.getName();
    m_array.resize(m_m*m_n);
    for (auto i = 0; i < m_m*m_n; i++) {
        m_array[i] = t_matrix(i);
    }
    return *this;
  }
  template<typename T>
  bool Matrix<T>::operator==(const Matrix<T>& t_matrix) const
  {
    if (m_n != t_matrix.getNumRows() || m_m != t_matrix.getNumCols()) {
      return false;
    }
    for (auto i = 0; i < m_m*m_n; i++) {
        if (t_matrix(i) != m_array[i]) {
          return false;
        }
    }
    return true;
  }
  template<typename T>
  bool Matrix<T>::operator!=(const Matrix<T>& t_matrix) const
  {
    if (m_n != t_matrix.getNumRows() || m_m != t_matrix.getNumCols()) {
      return true;
    }
    for (auto i = 0; i < m_m*m_n; i++) {
        if (t_matrix(i) != m_array[i]) {
          return true;
        }
    }
    return false;
  }
  template<typename T>
  Matrix<T> Matrix<T>::operator-() const
  {
    std::vector<T> mat(m_m*m_n);
    for (auto i = 0; i < m_m*m_n; i++) {
      mat[i] = -1*m_array[i];
    }
    return Matrix<T>(m_name,m_m,m_n,mat);
  }
  //  Matrix algebra
  template<typename T>
  Matrix<T> Matrix<T>::operator+(const Matrix<T>& t_matrix) const
  {
    if(m_n != t_matrix.getNumRows()) {
      Matrix<T> copy(*this);
      copy.setInfo(MATRIX_ADD_INCOMPATIBLE_ROWS(m_n, t_matrix.getNumRows(),
                                           m_name, t_matrix.getName()));
      return copy;
    }
    if (m_m != t_matrix.getNumCols()) {
      Matrix<T> copy(*this);
      copy.setInfo(MATRIX_ADD_INCOMPATIBLE_ROWS(m_n, t_matrix.getNumRows(),
                                           m_name, t_matrix.getName()));
      return copy;
    }
    std::string name = "(" + m_name + " + " + t_matrix.getName() + ")";
    Matrix<T> l(name, m_m, m_n, 0.0);
    for (auto i = 0; i < m_m*m_n; i++) {
        l(i) = m_array[i] + t_matrix(i);
    }
    return l;
  }
  template<typename T>
  Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& t_matrix)
  {
    if(m_n != t_matrix.getNumRows() || m_m != t_matrix.getNumCols()) {
      std::cout << "Matrices incompatible!" << std::endl;
      return *this;
    }
    std::string name = "(" + m_name + " + " + t_matrix.getName() + ")";
    setName(name);
    for (auto i = 0; i < m_m*m_n; i++) {
        m_array[i] += t_matrix(i);
    }
    return *this;
  }
  template<typename T>
  Matrix<T> Matrix<T>::operator-(const Matrix<T>& t_matrix) const
  {
    if(m_n != t_matrix.getNumRows() || m_m != t_matrix.getNumCols()) {
      std::cout << "Matrices incompatible!" << std::endl;
      return *this;
    }
    std::string name = "(" + m_name + " - " + t_matrix.getName() + ")";
    Matrix<T> l(name,m_m,m_n,0.0);
    for (auto i = 0; i < m_m*m_n; i++) {
        l(i) = m_array[i] - t_matrix(i);
    }
    return l;
  }
  template<typename T>
  Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& t_matrix)
  {
    if(m_n != t_matrix.getNumRows() || m_m != t_matrix.getNumCols()) {
      std::cout << "Matrices incompatible!" << std::endl;
      return *this;
    }
    std::string name = "(" + m_name + " - " + t_matrix.getName() + ")";
    setName(name);
    for (auto i = 0; i < m_m*m_n; i++) {
        m_array[i] -= t_matrix(i);
    }
    return *this;
  }
  template<typename T>
  Matrix<T> Matrix<T>::operator*(const Matrix<T>& t_matrix) const
  {
    if(m_n != t_matrix.getNumRows()) {
      std::cout << "Matrices " + m_name + " and " + t_matrix.getName()
          + " incompatible!" << std::endl;
      return *this;
    }
    return DGEMM(*this,t_matrix);
  }
  template<typename T>
  Matrix<T> Matrix<T>::brutem_mul(const Matrix<T>& t_matrix) const
  {
    if(m_n != t_matrix.getNumRows()) {
      std::cout << "Matrices incompatible!" << std::endl;
      return *this;
    }
    std::string name = "(" + m_name + " * " + t_matrix.getName() + ")";
    Matrix<T> l(name,m_m,t_matrix.getNumCols(),0.0);
    for (auto i = 0; i < m_m; i++) {
      for (auto j = 0; j < m_n; j++) {
        for (auto k = 0; k < m_m; k++) {
          l(i,j) += this->m_array[i*m_n + k] * t_matrix(k,j);
        }
      }
    }
    return l;
  }
  template<typename T>
  Matrix<T>& Matrix<T>::operator*=(const Matrix<T>& t_matrix)
  {
    Matrix<T> l = (*this) * t_matrix;
    (*this) = l;
    return *this;
  }
  //  Scalar algebra
  template<typename T>
  Matrix<T> Matrix<T>::operator+(const T& t_s) const
  {
    std::string name = "(" + m_name + " + " + std::to_string(t_s) + "I)";
    Matrix<T> l(name,m_m,m_n,0.0);
    for (auto i = 0; i < m_m*m_n; i++)
    {
        l(i) = m_array[i] + t_s;
    }
    return l;
  }
  template<typename T>
  Matrix<T> Matrix<T>::operator-(const T& t_s) const
  {
    std::string name = "(" + m_name + " - " + std::to_string(t_s) + "I)";
    Matrix<T> l(name,m_m,m_n,0.0);
    for (auto i = 0; i < m_m*m_n; i++)
    {
        l(i) = m_array[i] - t_s;
    }
    return l;
  }
  template<typename T>
  Matrix<T> Matrix<T>::operator*(const T& t_s) const
  {
    std::string name = "(" + m_name + " * " + std::to_string(t_s) + ")";
    Matrix<T> l(name,m_m,m_n,0.0);
    for (auto i = 0; i < m_m*m_n; i++)
    {
        l(i) = m_array[i] * t_s;
    }
    return l;
  }
  template<typename T>
  Matrix<T> Matrix<T>::operator/(const T& t_s) const
  {
    std::string name = "(" + m_name + " / " + std::to_string(t_s) + ")";
    Matrix<T> l(name,m_m,m_n,0.0);
    if (t_s == 0) {
      return l;
    }
    for (auto i = 0; i < m_m*m_n; i++) {
        l(i) = m_array[i] + t_s;
    }
    return l;
  }

  template<typename T>
  Matrix<T>& Matrix<T>::operator+=(const T& t_s)
  {
    std::string name = "(" + m_name + " + " + std::to_string(t_s) + "I)";
    setName(name);
    for (auto i = 0; i < m_m*m_n; i++) {
        m_array[i] += t_s;
    }
    return *this;
  }
  template<typename T>
  Matrix<T>& Matrix<T>::operator-=(const T& t_s)
  {
    std::string name = "(" + m_name + " - " + std::to_string(t_s) + "I)";
    setName(name);
    for (auto i = 0; i < m_m*m_n; i++) {
        m_array[i] -= t_s;
    }
    return *this;
  }
  template<typename T>
  Matrix<T>& Matrix<T>::operator*=(const T& t_s)
  {
    std::string name = "(" + m_name + " * " + std::to_string(t_s) + ")";
    setName(name);
    for (auto i = 0; i < m_m*m_n; i++) {
        m_array[i] *= t_s;
    }
    return *this;
  }
  template<typename T>
  Matrix<T>& Matrix<T>::operator/=(const T& t_s)
  {
    if (t_s == 0) {
      return *this;
    }
    std::string name = "(" + m_name + " / " + std::to_string(t_s) + ")";
    setName(name);
    for (auto i = 0; i < m_m*m_n; i++) {
        m_array[i] += t_s;
    }
    return *this;
  }
  //  Multiplying a vector
  template<typename T>
  Vector<T> Matrix<T>::operator*(const Vector<T>& t_v)
  {
    std::vector<T> vec(m_n,0.0);
    Vector<T> v2(vec);
    if (m_n != t_v.getDim()) {
      return v2;
    }
    for (auto i = 0; i < m_m; i++) {
      T temp = 0.0;
      for (auto j = 0; j < m_n; j++) {
        temp += this->m_array[i*m_n + j] * t_v(j);
      }
      v2(i) = temp;
    }
    return v2;
  }
  //  Access Matrix<T>::operators
  template<typename T>
  T& Matrix<T>::operator()(const size_t& t_i, const size_t& t_j)
  {
    return m_array[t_i*m_n + t_j];
  }
  template<typename T>
  const T& Matrix<T>::operator()(const size_t& t_i, const size_t& t_j) const
  {
    return m_array[t_i*m_n + t_j];
  }
  template<typename T>
  T& Matrix<T>::operator()(const size_t& t_i)
  {
    return m_array[t_i];
  }
  template<typename T>
  const T& Matrix<T>::operator()(const size_t& t_i) const
  {
    return m_array[t_i];
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Various methods which include,
  //    print()
  //    summary()
  //    tranpose()
  //    transpose_inplace()
  //    trace()
  //----------------------------------------------------------------------------
  template<typename T>
  void Matrix<T>::print()
  {
    std::cout << "(" << m_m << " x " << m_n << ") Matrix";
    if (m_name != " ") {
      std::cout << ": '" << m_name << "'";
    }
    std:: cout << "\n[ ";
    for (auto i = 0; i < m_n; i++) {
      for (auto j = 0; j < m_m; j++) {
        std::cout << this->m_array[i*m_m + j];
        if (j < m_m-1)
        std::cout << "   ";
      }
      if (i < m_n-1)
        std::cout << "\n  ";
    }
    std::cout << " ]" << std::endl;
  }
  template<typename T>
  const std::string Matrix<T>::summary()
  {
    std::stringstream s;
    s.str("");
    s.clear();
    std::string sum = "dim: (" + std::to_string(m_m)
                    + "x" + std::to_string(m_n)
                    + "), type: "
                    + type_name<decltype(m_array[0])>();
    if (m_name != " ") {
      sum +=  ", name: '" + m_name + "'";
    }
    if (m_array.size() == 0) {
      sum += "\n[  empty  ]";
      return sum;
    }
    sum += "\n[ ";
    if (m_m < 10) {
      for (auto i = 0; i < m_m; i++) {
        if (m_n < 10) {
          if (m_array[i*m_n] >= 0.0) {
            sum += " ";
          }
          for (auto j = 0; j < m_n; j++) {
            sum += scientific_not(this->m_array[i*m_n + j],3);
            if (j < m_n-1) {
              if (m_array[i*m_n + j+1] >= 0.0) {
                sum += "   ";
              }
              else {
                sum += "  ";
              }
            }
          }
        }
        else {
          if (m_array[i*m_n] >= 0.0) {
            sum += " ";
          }
          sum += scientific_not(this->m_array[i*m_n + 0],3);
          if (m_array[i*m_n + 1] >= 0.0) {
            sum += "   ";
          }
          else {
            sum += "  ";
          }
          sum += scientific_not(this->m_array[i*m_n + 1],3);
          if (m_array[i*m_n + 2] >= 0.0) {
            sum += "   ";
          }
          else {
            sum += "  ";
          }
          sum += scientific_not(this->m_array[i*m_n + 2],3);
          sum += "   ";
          sum += "...   ...   ...   ";
          if (m_array[i*m_n + m_n-3] >= 0.0) {
            sum += " ";
          }
          sum += scientific_not(this->m_array[i*m_n + m_n-3],3);
          if (m_array[i*m_n + m_n-2] >= 0.0) {
            sum += "   ";
          }
          else {
            sum += "  ";
          }
          sum += scientific_not(this->m_array[i*m_n + m_n-2],3);
          if (m_array[i*m_n + m_n-1] >= 0.0) {
            sum += "   ";
          }
          else {
            sum += "  ";
          }
          sum += scientific_not(this->m_array[i*m_n + m_n-1],3);
        }
        if (i < m_m-1) {
          sum += "\n  ";
        }
      }
    }
    else {
      for (auto i = 0; i < 3; i++) {
        if (m_n < 10) {
          if (m_array[i*m_n] >= 0.0) {
            sum += " ";
          }
          for (auto j = 0; j < m_n; j++) {
            sum += scientific_not(this->m_array[i*m_n + j],3);
            if (j < m_n-1) {
              if (m_array[i*m_n + j+1] >= 0.0) {
                sum += "   ";
              }
              else {
                sum += "  ";
              }
            }
          }
        }
        else {
          if (m_array[i*m_n] >= 0.0) {
            sum += " ";
          }
          sum += scientific_not(this->m_array[i*m_n + 0],3);
          if (m_array[i*m_n + 1] >= 0.0) {
            sum += "   ";
          }
          else {
            sum += "  ";
          }
          sum += scientific_not(this->m_array[i*m_n + 1],3);
          if (m_array[i*m_n + 2] >= 0.0) {
            sum += "   ";
          }
          else {
            sum += "  ";
          }
          sum += scientific_not(this->m_array[i*m_n + 2],3);
          sum += "   ";
          sum += "...   ...   ...   ";
          if (m_array[i*m_n + m_n-3] >= 0.0) {
            sum += " ";
          }
          sum += scientific_not(this->m_array[i*m_n + m_n-3],3);
          if (m_array[i*m_n + m_n-2] >= 0.0) {
            sum += "   ";
          }
          else {
            sum += "  ";
          }
          sum += scientific_not(this->m_array[i*m_n + m_n-2],3);
          if (m_array[i*m_n + m_n-1] >= 0.0) {
            sum += "   ";
          }
          else {
            sum += "  ";
          }
          sum += scientific_not(this->m_array[i*m_n + m_n-1],3);
        }
        if (i < m_m-1) {
          sum += "\n  ";
        }
      }
      sum += "    ...";
      if (m_n < 10) {
        for (auto j = 0; j < m_n-1; j++) {
          sum += "         ...";
        }
        sum += "\n  ";
      }
      else {
        sum += "         ...         ...      ...   ...   ...       ...         ...         ...\n  ";
      }
      for (auto i = m_m-3; i < m_m; i++) {
        if (m_n < 10) {
          if (m_array[i*m_n] >= 0.0) {
            sum += " ";
          }
          for (auto j = 0; j < m_n; j++) {
            sum += scientific_not(this->m_array[i*m_n + j],3);
            if (j < m_n-1) {
              if (m_array[i*m_n + j+1] >= 0.0) {
                sum += "   ";
              }
              else {
                sum += "  ";
              }
            }
          }
        }
        else {
          if (m_array[i*m_n] > 0.0) {
            sum += " ";
          }
          sum += scientific_not(this->m_array[i*m_n + 0],3);
          if (m_array[i*m_n + 1] >= 0.0) {
            sum += "   ";
          }
          else {
            sum += "  ";
          }
          sum += scientific_not(this->m_array[i*m_n + 1],3);
          if (m_array[i*m_n + 2] >= 0.0) {
            sum += "   ";
          }
          else {
            sum += "  ";
          }
          sum += scientific_not(this->m_array[i*m_n + 2],3);
          sum += "   ";
          sum += "...   ...   ...   ";
          if (m_array[i*m_n + m_n-3] >= 0.0) {
            sum += " ";
          }
          sum += scientific_not(this->m_array[i*m_n + m_n-3],3);
          if (m_array[i*m_n + m_n-2] >= 0.0) {
            sum += "   ";
          }
          else {
            sum += "  ";
          }
          sum += scientific_not(this->m_array[i*m_n + m_n-2],3);
          if (m_array[i*m_n + m_n-1] >= 0.0) {
            sum += "   ";
          }
          else {
            sum += "  ";
          }
          sum += scientific_not(this->m_array[i*m_n + m_n-1],3);
        }
        if (i < m_m-1) {
          sum += "\n  ";
        }
      }
    }
    sum += "  ]";
    return sum;
  }
  template<typename T>
  Matrix<T> Matrix<T>::transpose() const
  {
    std::vector<T> new_array(m_m*m_n);
    for (auto i = 0 ; i < m_m; i++)
    {
      for (auto j = 0; j < m_n; j++)
      {
        new_array[j*m_m + i] = m_array[i*m_n + j];
      }
    }
    std::string name = "(" + m_name + ")^T";
    int m = m_n;
    int n = m_m;
    return Matrix<T>(name,m,n,new_array);
  }
  template<typename T>
  void Matrix<T>::transpose_inplace(bool inplace)
  {
    std::vector<T> new_array(m_m*m_n);
    for (auto i = 0 ; i < m_m; i++)
    {
      for (auto j = 0; j < m_n; j++)
      {
        new_array[j*m_m + i] = m_array[i*m_n + j];
      }
    }
    std::string name = "(" + m_name + ")^T";
    m_name = name;
    int m = m_n;
    int n = m_m;
    m_n = n;
    m_m = m;
    m_array = new_array;
  }
  template<typename T>
  T Matrix<T>::trace()
  {
    if (m_m != m_n) {
      std::cout << "Matrix is not square, trace is undefined!" << std::endl;
      return 0;
    }
    T result = 0;
    for (auto i = 0; i < m_m; i++) {
      result += m_array[i*m_n + i];
    }
    return result;
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Various instantiators
  //----------------------------------------------------------------------------
  Matrix<double> identity_d(size_t t_m)
  {
    std::string name = "I_{" + std::to_string(t_m) + "x"
                     + std::to_string(t_m) + "}";
    Matrix<double> matrix(name, t_m, t_m, 0.0);

    for (auto i = 0; i < t_m; i++) {
      matrix(i,i) = 1.0;
    }
    return matrix;
  }
  Matrix<double> zeros_d(size_t t_m)
  {
    Matrix<double> z(t_m, t_m, 0.0);
    return z;
  }
  Matrix<double> zeros_d(size_t t_m, size_t t_n)
  {
    Matrix<double> z(t_m, t_n, 0.0);
    return z;
  }
  Matrix<double> ones_d(size_t t_m)
  {
    Matrix<double> o(t_m, t_m, 1.0);
    return o;
  }
  Matrix<double> ones_d(size_t t_m, size_t t_n)
  {
    Matrix<double> o(t_m, t_n, 1.0);
    return o;
  }
  //----------------------------------------------------------------------------
  //  Permutation matrix
  //----------------------------------------------------------------------------
  Matrix<double> permutationMatrix_d(const size_t& t_m,
                              const std::vector<size_t> t_pivot)
  {
    //  generate a permutation matrix from a set of pivot indices
    std::vector<double> swaps(t_m*t_m, 0);
    std::vector<size_t> p(t_m,0);
    for (int i = 0; i < t_m; i++) {
      p[i] = i;
    }
    for (int i = 0; i < t_m; i++) {
      int temp;
      temp = p[t_pivot[i]-1];
      p[t_pivot[i]-1] = p[i];
      p[i] = temp;
    }
    for(size_t i = 0; i < t_m; i++) {
      swaps[p[i]*t_m + i] = 1;
    }
    std::string name = "(" + std::to_string(t_m) = "x"
                        + std::to_string(t_m) + ") perm";
    return Matrix<double>(name,t_m,t_m,swaps);
  }
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Overload of the ostream operator<<
  //  This has the same functionality as print, except it can be used in
  //  a std::cout statement.
  //----------------------------------------------------------------------------
  template<typename T>
  std::ostream& operator<<(std::ostream& t_os, const Matrix<T>& t_matrix)
  {
    t_os << "(" << t_matrix.getNumRows() << " x ";
    t_os << t_matrix.getNumCols() << ") Matrix";
    if (t_matrix.getName() != " ") {
      t_os << ": '" << t_matrix.getName() << "'";
    }
    t_os << "\n[ ";
    for (auto i = 0; i < t_matrix.getNumRows(); i++) {
      t_os << "[ ";
      for (auto j = 0; j < t_matrix.getNumCols(); j++) {
        t_os << t_matrix(i,j) << " ";
      }
      t_os << "]";
      if (i < t_matrix.getNumRows()-1) {
        t_os << "\n  ";
      }
    }
    t_os << " ]" << std::endl;
    return t_os;
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Level 2 BLAS methods
  //  Matrix-vector multiplication
  //----------------------------------------------------------------------------
  //  DGEMV - generic matrix vector multiplication
  //  Arguments:  alpha - double
  //              A     - (m x n) matrix
  //              x     - (n)-dim vector
  //
  //  Returns:    alpha * A * x
  //----------------------------------------------------------------------------
  Vector<double> DGEMV(double& t_alpha, Matrix<double>& t_A,
                       Vector<double>& t_x)
  {
    //  container for the vector y
    std::vector<double> y(t_x.getDim());
    cblas_dgemv(CblasRowMajor,   //  row major order
                CblasNoTrans,    //  don't transpose A
                t_A.getNumRows(),//  number of rows of A
                t_A.getNumCols(),//  number of columns of A
                t_alpha,         //  scaling factor for A*x
                t_A.data(),      //  pointer to A t_array
                t_A.getNumRows(),//  number of rows in A
                t_x.data(),      //  pointer to x vector
                1,               //  increment for elements of x
                1,               //  scaling factor for y
                y.data(),        //  pointer to y vector
                1);              //  increment for elements of y
    //  generate a new name for the product
    std::string name;
    if (t_alpha != 0) {
      name += std::to_string(t_alpha) + " * ";
    }
    name += t_A.getName() + " * " + t_x.getName();
    Vector<double> result(name,y);
    return result;
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DGEMV with (beta,y) - generic matrix vector multiplication
  //  Arguments:  alpha - double
  //              A     - (m x n) matrix
  //              x     - (n)-dim vector
  //              beta  - double
  //              y     - (n)-dim vector
  //
  //  Returns:    alpha * A * x + beta * y
  //  This method overwrites the vector y!
  //----------------------------------------------------------------------------
  void DGEMV(double& t_alpha, Matrix<double>& t_A,
             Vector<double>& t_x, double& t_beta, Vector<double>& t_y)
  {
    cblas_dgemv(CblasRowMajor,   //  row major order
                CblasNoTrans,    //  don't transpose A
                t_A.getNumRows(),//  number of rows of A
                t_A.getNumCols(),//  number of columns of A
                t_alpha,         //  scaling factor for A*x
                t_A.data(),      //  pointer to A t_array
                t_A.getNumRows(),//  number of rows in A
                t_x.data(),      //  pointer to x vector
                1,               //  increment for elements of x
                t_beta,          //  scaling factor for y
                t_y.data(),      //  pointer to y data
                1);              //  increment for the elements of y
    //  generate a new name for the product
    std::string name;
    if (t_alpha != 0) {
      name += std::to_string(t_alpha) + " * ";
    }
    name += t_A.getName() + " * " + t_x.getName();
    if (t_beta != 0) {
      name += " + " + std::to_string(t_beta) + " * " + t_y.getName();
    }
    t_y.setName(name);
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DGER - generic outer-product of two vectors
  //  Arguments:  alpha - double
  //              x     - (m)-dim vector
  //              y     - (n)-dim vector
  //
  //  Returns:    alpha * x * y^T ( (m x n) - matrix )
  //----------------------------------------------------------------------------
  Matrix<double> DGER(double& t_alpha, Vector<double>& t_x, Vector<double>& t_y)
  {
    std::vector<double> A(t_x.getDim() * t_y.getDim());
    cblas_dger(CblasRowMajor,  //  row major order
               t_x.getDim(),   //  number of rows of A
               t_y.getDim(),   //  number of columns of A
               t_alpha,        //  scaling factor for A*x
               t_x.data(),     //  pointer to x vector
               1,              //  increment for elements of x
               t_y.data(),     //  pointer to y data
               1,              //  increment for the elements of y
               A.data(),       //  pointer to the matrix A
               t_y.getDim());  //  leading dimension of A
    //  generate a new name for the product
    std::string name;
    if (t_alpha != 0 && t_x.getName() != " " && t_y.getName() != " ") {
      name += std::to_string(t_alpha) + " * " + t_x.getName()
              + " * " + t_y.getName() + "^T";
    }
    return Matrix<double>(name,t_x.getDim(),t_y.getDim(),A);
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DGER (in place) - generic outer-product of two vectors
  //  Arguments:  alpha - double
  //              x     - (m)-dim vector
  //              y     - (n)-dim vector
  //              A     - (m x n) Matrix
  //
  //  Returns:    void
  //----------------------------------------------------------------------------
  void DGER(double& t_alpha, Vector<double>& t_x, Vector<double>& t_y,
                      Matrix<double>& t_matrix)
  {
    if (t_x.getDim() != t_matrix.getNumRows() ||
        t_y.getDim() != t_matrix.getNumCols())
    {
      std::cout << "Vectors and Matrix are incompatible!" << std::endl;
      return;
    }
    cblas_dger(CblasRowMajor,//  row major order
               t_x.getDim(),   //  number of rows of A
               t_y.getDim(),   //  number of columns of A
               t_alpha,        //  scaling factor for A*x
               t_x.data(),     //  pointer to x vector
               1,            //  increment for elements of x
               t_y.data(),     //  pointer to y data
               1,            //  increment for the elements of y
               t_matrix.data(),     //  pointer to the matrix A
               t_y.getDim());  //  leading dimension of A
    //  generate a new name for the product
    std::string name;
    if (t_alpha != 0 && t_x.getName() != " " && t_y.getName() != " ") {
      name += std::to_string(t_alpha) + " * " + t_x.getName()
              + " * " + t_y.getName() + "^T";
    }
    if (t_matrix.getName() != " ") {
      name += " + " + t_matrix.getName();
    }
    t_matrix.setName(name);
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Level 3 BLAS methods
  //  Matrix-matrix multiplication
  //----------------------------------------------------------------------------
  //  DGEMM - generic matrix-matrix multiplication
  //  Arguments:  alpha - double
  //              A     - (m x k)-matrix
  //              B     - (k x n)-matrix
  //
  //  Returns:    alpha * m * n  ( (m x n)-matrix )
  //----------------------------------------------------------------------------
  Matrix<double> DGEMM(const double& t_alpha, const Matrix<double>& t_A,
                       const Matrix<double>& t_B)
  {
    if (t_A.getNumCols() != t_B.getNumRows())
    {
      std::cout << "Matrices are incompatible!" << std::endl;
      return Matrix<double>("zeros",1,0.0);
    }
    Matrix<double> A_copy(t_A);
    Matrix<double> B_copy(t_B);
    std::vector<double> C(t_A.getNumRows() * t_B.getNumCols());
    cblas_dgemm(CblasRowMajor,    //  row major order
                CblasNoTrans,     //  don't tranpose A
                CblasNoTrans,     //  don't transpose B
                t_A.getNumRows(), //  number of rows of A
                t_B.getNumCols(), //  number of cols of B
                t_A.getNumCols(), //  number of cols of A
                t_alpha,          //  scaling factor for A*B
                A_copy.data(),    //  pointer to elements of A
                t_A.getNumCols(), //  leading dimension of A
                B_copy.data(),    //  pointer to elements of B
                t_B.getNumCols(), //  leading dimension of B
                1,                //  coefficient beta=1
                C.data(),         //  pointer to elements of C
                t_B.getNumCols());//  leading dimension of C
    //  generate a new name for the product
    std::string name;
    if (t_alpha != 0 && t_A.getName() != " " && t_B.getName() != " ") {
      name += std::to_string(t_alpha) + " * " + t_A.getName()
              + " * " + t_B.getName();
    }
    return Matrix<double>(name,t_A.getNumRows(),t_B.getNumCols(),C);
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DGEMM - generic matrix-matrix multiplication
  //  Arguments:  alpha - double
  //              A     - (m x k)-matrix
  //              B     - (k x n)-matrix
  //              beta  - double
  //              C     - (m x n)-matrix
  //
  //  Returns:    alpha * m * n  ( (m x n)-matrix )
  //----------------------------------------------------------------------------
  void DGEMM(const double& t_alpha, const Matrix<double>& t_A,
             const Matrix<double>& t_B, const double& t_beta,
             Matrix<double>& t_C)
  {
    if (t_A.getNumCols() != t_B.getNumRows() ||
        t_A.getNumRows() != t_C.getNumRows() ||
        t_B.getNumCols() != t_C.getNumCols())
    {
      std::cout << "Matrices are incompatible!" << std::endl;
      return;
    }
    Matrix<double> A_copy(t_A);
    Matrix<double> B_copy(t_B);
    cblas_dgemm(CblasRowMajor,    //  row major order
                CblasNoTrans,     //  don't tranpose A
                CblasNoTrans,     //  don't transpose B
                t_A.getNumRows(), //  number of rows of A
                t_B.getNumCols(), //  number of cols of B
                t_A.getNumCols(), //  number of cols of A
                t_alpha,          //  scaling factor for A*B
                A_copy.data(),    //  pointer to elements of A
                t_A.getNumCols(), //  leading dimension of A
                B_copy.data(),    //  pointer to elements of B
                t_B.getNumCols(), //  leading dimension of B
                t_beta,           //  coefficient beta=1
                t_C.data(),       //  pointer to elements of C
                t_B.getNumCols());//  leading dimension of C
    //  generate a new name for the product
    std::string name;
    if (t_alpha != 0 && t_A.getName() != " " && t_B.getName() != " ") {
      name += std::to_string(t_alpha) + " * " + t_A.getName()
              + " * " + t_B.getName();
    }
    t_C.setName(name);
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DGEMM (no alpha) - generic matrix-matrix multiplication
  //  Arguments:  A     - (m x k)-matrix
  //              B     - (k x n)-matrix
  //
  //  Returns:    m * n  ( (m x n)-matrix )
  //----------------------------------------------------------------------------
  Matrix<double> DGEMM(const Matrix<double>& t_A,
                       const Matrix<double>& t_B)
  {
    if (t_A.getNumCols() != t_B.getNumRows()) {
      std::cout << "Matrices are incompatible!" << std::endl;
      return Matrix<double>("zeros",1,0.0);
    }
    Matrix<double> A_copy(t_A);
    Matrix<double> B_copy(t_B);
    std::vector<double> C(t_A.getNumRows() * t_B.getNumCols());
    cblas_dgemm(CblasRowMajor,    //  row major order
                CblasNoTrans,     //  don't tranpose A
                CblasNoTrans,     //  don't transpose B
                t_A.getNumRows(), //  number of rows of A
                t_B.getNumCols(), //  number of cols of B
                t_A.getNumCols(), //  number of cols of A
                1.0,              //  scaling factor for A*B
                A_copy.data(),    //  pointer to elements of A
                t_A.getNumCols(), //  leading dimension of A
                B_copy.data(),    //  pointer to elements of B
                t_B.getNumCols(), //  leading dimension of B
                1.0,              //  coefficient beta=1
                C.data(),         //  pointer to elements of C
                t_B.getNumCols());//  leading dimension of C
    //  generate a new name for the product
    std::string name;
    if (t_A.getName() != " " && t_B.getName() != " ") {
      name += t_A.getName() + " * " + t_B.getName();
    }
    return Matrix<double>(name,t_A.getNumRows(),t_B.getNumCols(),C);
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  LAPACK methods
  //----------------------------------------------------------------------------
  //  DGELS - solve the linear system min_C|B - AC|
  //  Arguments:  A     - (m x k)-matrix
  //              B     - (m x n)-matrix
  //
  //  Returns:    Cm_min  ( (k x n)-matrix )
  //----------------------------------------------------------------------------
  Matrix<double> DGELS(const Matrix<double>& t_A, const Matrix<double>& t_B)
  {
    if (t_A.getNumRows() != t_B.getNumRows()) {
      std::cout << "Matrices are incompatible!" << std::endl;
      return Matrix<double>("zeros",1,0.0);
    }
    Matrix<double> QR(t_A);
    std::vector<double> c = t_B.getArray();
    int info;
    info = LAPACKE_dgels(LAPACK_ROW_MAJOR,  //  row major layout
                         'N',               //  don't transpose A
                         t_A.getNumRows(),  //  number of rows of A
                         t_A.getNumCols(),  //  number of columns of A
                         t_B.getNumCols(),  //  number of columns of B
                         QR.data(),         //  pointer to elements of A
                         t_A.getNumCols(),  //  leading dimension of A
                         c.data(),          //  pointer to elements of B
                         t_B.getNumCols()); //  leading dimension of B
    //  Cut the result according to (A.getNumCols() x B.getNumCols())
    c.resize(t_A.getNumCols() * t_B.getNumCols());
    std::string name;
    if (t_A.getName() != " " && t_B.getName() != " ") {
      name += "min_C||" + t_A.getName() + "*C - " + t_B.getName() + "|";
    }
    else {
      name  = " ";
    }
    Matrix<double> C(name,t_A.getNumCols(),t_B.getNumCols(),c);
    if( info > 0 ) {
      std::string s = "The diagonal element " + std::to_string(info)
                    + " of the triangular factor of A is zero, so that A "
                    + "does not have full rank;the least squares solution"
                    + " could not be computed.";
      C.setFlag(-1);
      C.setInfo(s);
    }
    else if (info < 0) {
      C.setFlag(-1);
      C.setInfo(std::to_string(info)+"th argument has illegal value");
    }
    return C;
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DGELS - solve the linear system min_u|v - Au|
  //  Arguments:  A     - (m x k)-matrix
  //              v     - (m)-dim vector
  //
  //  Returns:    um_min  ( (k)-dim vector )
  //----------------------------------------------------------------------------
  Vector<double> DGELS(const Matrix<double>& t_A, const Vector<double>& t_v)
  {
    if (t_A.getNumRows() != t_v.getDim())
    {
      std::cout << "Matrices and vector are incompatible!" << std::endl;
      return Vector<double>("zeros",1,0.0);
    }
    Matrix<double> QR(t_A);
    //  ------------------------------------------------------------------------
    //  WARNING: while A can be a general m x n matrix, the vector u
    //  must have at least as many rows as the number of columns of A,
    //  otherwise we get a segfault from LAPACK.  To prevent this, u
    //  is t_initialized with _dim = A.getNumCols() with zeros.  Then, the
    //  first v.getDim() elements are filled with v's elements.
    //  (N. Carrara - 6/30/2020)
    //  Error fixed with commit - bc555791c2ebc9aadebe44e5b74fbc367b7e2123.
    //--------------------------------------------------------------------------
    //  WARNING: The previous correction was actually wrong.  The vector u
    //  must be t_initialized with std::max(A.getNumCols(),A.getNumRows()).
    //  (N. Carrara - 6/30/2020)
    //  Error fixed with commit - ab0767e605414d121ff6c8d844afc28e9867ea0f.
    //--------------------------------------------------------------------------
    std::vector<double> u(std::max(t_A.getNumCols(),t_A.getNumRows()),0.0);
    for (auto i = 0; i < t_v.getDim(); i++) {
      u[i] = t_v(i);
    }
    //std::cout << "\nC++ Location of v: " << &v;
    int info;
    info = LAPACKE_dgels(LAPACK_ROW_MAJOR,//  row major layout
                         'N',             //  don't transpose A
                         t_A.getNumRows(),  //  number of rows of A
                         t_A.getNumCols(),  //  number of rows of A
                         1,               //  dimension of v
                         QR.data(),       //  pointer to elements of A
                         t_A.getNumCols(),  //  leading dimension of A
                         u.data(),        //  pointer to elements of v
                         1);              //  leading dimension of v
    //  Cut the result according to (A.getNumCols())
    u.resize(t_A.getNumCols());
    std::string name;
    if (t_A.getName() != " " && t_v.getName() != " ") {
      name += "min_u||" + t_A.getName() + "*u - " + t_v.getName() + "|";
    }
    else {
      name  = " ";
    }
    Vector<double> U(name,u);
    if( info > 0 ) {
     std::string s = "The diagonal element " + std::to_string(info)
                   + " of the triangular factor of A is zero, so that A "
                   + "does not have full rank;the least squares solution"
                   + " could not be computed.";
     U.setFlag(-1);
     U.setInfo(s);
    }
    else if (info < 0) {
     U.setFlag(-1);
     U.setInfo(std::to_string(info)+"th argument has illegal value");
    }
    return U;
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DGELSY - solve the linear system min_C|B - AC| using
  //           complete orthogonal factorization.
  //  Arguments:  A     - (m x k)-matrix
  //              B     - (m x n)-matrix
  //
  //  Returns:    Cm_min  ( (k x n)-matrix )
  //----------------------------------------------------------------------------
  Matrix<double> DGELSY(const Matrix<double>& t_A, const Matrix<double>& t_B)
  {
    if (t_A.getNumRows() != t_B.getNumRows()) {
      std::cout << "Matrices are incompatible!" << std::endl;
      return Matrix<double>("zeros",1,0.0);
    }
    Matrix<double> QR(t_A);
    std::vector<double> c = t_B.getArray();
    int info;
    int jpvt[t_A.getNumCols()];
    double rcond;
    int rank;
    info = LAPACKE_dgelsy(LAPACK_ROW_MAJOR,  //  row major layout
                          t_A.getNumRows(),  //  number of rows of A
                          t_A.getNumCols(),  //  number of columns of A
                          t_B.getNumCols(),  //  number of columns of B
                          QR.data(),         //  pointer to elements of A
                          t_A.getNumCols(),  //  leading dimension of A
                          c.data(),          //  pointer to elements of B
                          t_B.getNumCols(),  //  leading dimension of B
                          jpvt,              //  workspace t_array
                          rcond,             //  used for rank of A
                          &rank);            //  the effective rank of A
    //  Cut the result according to (A.getNumCols() x B.getNumCols())
    c.resize(t_A.getNumCols() * t_B.getNumCols());
    std::string name;
    if (t_A.getName() != " " && t_B.getName() != " ") {
      name += "min_C||" + t_A.getName() + "*C - " + t_B.getName() + "|";
    }
    else {
      name  = " ";
    }
    Matrix<double> C(name,t_A.getNumCols(),t_B.getNumCols(),c);
    if( info > 0 ) {
      std::string s = "The diagonal element " + std::to_string(info)
                    + " of the triangular factor of A is zero, so that A "
                    + "does not have full rank;the least squares solution"
                    + " could not be computed.";
      C.setFlag(-1);
      C.setInfo(s);
    }
    else if (info < 0) {
      C.setFlag(-1);
      C.setInfo(std::to_string(info)+"th argument has illegal value");
    }
    return C;
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DGELSY - solve the linear system min_u|v - Au|
  //           using complete orthogonal factorization
  //  Arguments:  A     - (m x k)-matrix
  //              v     - (m)-dim vector
  //A.getNumCols()
  //  Returns:    um_min  ( (k)-dim vector )
  //----------------------------------------------------------------------------
  Vector<double> DGELSY(const Matrix<double>& t_A, const Vector<double>& t_v)
  {
    if (t_A.getNumRows() != t_v.getDim()) {
      std::cout << "Matrices and vector are incompatible!" << std::endl;
      return Vector<double>("zeros",1,0.0);
    }
    Matrix<double> QR(t_A);
    //  ------------------------------------------------------------------------
    //  WARNING: while A can be a general m x n matrix, the vector u
    //  must have at least as many rows as the number of columns of A,
    //  otherwise we get a segfault from LAPACK.  To prevent this, u
    //  is t_initialized with _dim = A.getNumCols() with zeros.  Then, the
    //  first v.getDim() elements are filled with v's elements.
    //  (N. Carrara - 6/30/2020)
    //  Error fixed with commit - bc555791c2ebc9aadebe44e5b74fbc367b7e2123.
    //--------------------------------------------------------------------------
    std::vector<double> u(std::max(t_A.getNumCols(),t_A.getNumRows()),0.0);
    for (auto i = 0; i < t_v.getDim(); i++)
    {
      u[i] = t_v(i);
    }
    int info;
    int jpvt[t_A.getNumCols()];
    double rcond;
    int rank;
    info = LAPACKE_dgelsy(LAPACK_ROW_MAJOR,  //  row major layout
                          t_A.getNumRows(),  //  number of rows of A
                          t_A.getNumCols(),  //  number of columns of A
                          1,                 //  dimension of v
                          QR.data(),         //  pointer to elements of A
                          t_A.getNumCols(),  //  leading dimension of A
                          u.data(),          //  pointer to elements of v
                          1,                 //  leading dimension of v
                          jpvt,              //  workspace t_array
                          rcond,             //  used for rank of A
                          &rank);            //  the effective rank of A
    if( info > 0 ) {
      std::cout << "The diagonal element " << info << " of the triangular "
                << "factor of A is zero, so that A does not have full rank;\n"
                << "the least squares solution could not be computed.\n";
    }
    //  Cut the result according to (A.getNumCols())
    u.resize(t_A.getNumCols());
    std::string name;
    if (t_A.getName() != " " && t_v.getName() != " ") {
      name += "min_u||" + t_A.getName() + "*u - " + t_v.getName() + "|";
    }
    else {
      name  = " ";
    }
    Vector<double> U(name,u);
    if( info > 0 ) {
     std::string s = "The diagonal element " + std::to_string(info)
                   + " of the triangular factor of A is zero, so that A "
                   + "does not have full rank;the least squares solution"
                   + " could not be computed.";
     U.setFlag(-1);
     U.setInfo(s);
    }
    else if (info < 0) {
     U.setFlag(-1);
     U.setInfo(std::to_string(info)+"th argument has illegal value");
    }
    return U;
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DGELSD - solve the linear system min_C|B - AC|
  //           using SVD and the divide and conquer method
  //  Arguments:  A     - (m x k)-matrix
  //              B     - (m x n)-matrix
  //
  //  Returns:    Cm_min  ( (k x n)-matrix )
  //----------------------------------------------------------------------------
  Matrix<double> DGELSD(const Matrix<double>& t_A, const Matrix<double>& t_B)
  {
    if (t_A.getNumRows() != t_B.getNumRows()) {
      std::cout << "Matrices are incompatible!" << std::endl;
      return Matrix<double>("zeros",1,0.0);
    }
    Matrix<double> SVD(t_A);
    std::vector<double> c = t_B.getArray();
    int info;
    std::vector<double> singular(std::min(t_A.getNumRows(),t_A.getNumCols()));
    double rcond;
    int rank;
    info = LAPACKE_dgelsd(LAPACK_ROW_MAJOR,  //  row major layout
                          t_A.getNumRows(),  //  number of rows of A
                          t_A.getNumCols(),  //  number of columns of A
                          t_B.getNumCols(),  //  number of columns of B
                          SVD.data(),        //  pointer to elements of A
                          t_A.getNumCols(),  //  leading dimension of A
                          c.data(),          //  pointer to elements of B
                          t_B.getNumCols(),  //  leading dimension of B
                          singular.data(),   //  pointer to the singular values
                          rcond,             //  condition number
                          &rank);            //  effective rank of A.
    //  Cut the result according to (A.getNumCols() x B.getNumCols())
    c.resize(t_A.getNumCols() * t_B.getNumCols());
    std::string name;
    if (t_A.getName() != " " && t_B.getName() != " ") {
      name += "min_C||" + t_A.getName() + "*C - " + t_B.getName() + "|";
    }
    else {
      name  = " ";
    }
    Matrix<double> C(name,t_A.getNumCols(),t_B.getNumCols(),c);
    if( info > 0 ) {
      std::string s = "The diagonal element " + std::to_string(info)
                    + " of the triangular factor of A is zero, so that A "
                    + "does not have full rank;the least squares solution"
                    + " could not be computed.";
      C.setFlag(-1);
      C.setInfo(s);
    }
    else if (info < 0) {
      C.setFlag(-1);
      C.setInfo(std::to_string(info)+"th argument has illegal value");
    }
    return C;
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DGELSD - solve the linear system min_u|v - Au|
  //           using SVD and the divide and conquer method
  //  Arguments:  A     - (m x k)-matrix
  //              v     - (m)-dim vector
  //
  //  Returns:    um_min  ( (k)-dim vector )
  //----------------------------------------------------------------------------
  Vector<double> DGELSD(const Matrix<double>& t_A, const Vector<double>& t_v)
  {
    if (t_A.getNumRows() != t_v.getDim()) {
      std::cout << "Matrices and vector are incompatible!" << std::endl;
      return Vector<double>("zeros",1,0.0);
    }
    Matrix<double> SVD(t_A);
    //  ------------------------------------------------------------------------
    //  WARNING: while A can be a general m x n matrix, the vector u
    //  must have at least as many rows as the number of columns of A,
    //  otherwise we get a segfault from LAPACK.  To prevent this, u
    //  is t_initialized with _dim = A.getNumCols() with zeros.  Then, the
    //  first v.getDim() elements are filled with v's elements.
    //  (N. Carrara - 6/30/2020)
    //  Error fixed with commit - bc555791c2ebc9aadebe44e5b74fbc367b7e2123.
    //--------------------------------------------------------------------------
    std::vector<double> u(std::max(t_A.getNumCols(),t_A.getNumRows()),0.0);
    for (auto i = 0; i < t_v.getDim(); i++) {
      u[i] = t_v(i);
    }
    int info;
    std::vector<double> singular(std::min(t_A.getNumRows(),t_A.getNumCols()));
    double rcond;
    int rank;
    info = LAPACKE_dgelsd(LAPACK_ROW_MAJOR,  //  row major layout
                          t_A.getNumRows(),  //  number of rows of A
                          t_A.getNumCols(),  //  number of columns of A
                          1,                 //  dimension of v
                          SVD.data(),        //  pointer to elements of A
                          t_A.getNumCols(),  //  leading dimension of A
                          u.data(),          //  pointer to elements of v
                          1,                 //  leading dimension of v
                          singular.data(),   //  pointer to the singular values
                          rcond,             //  condition number
                          &rank);            //  effective rank of A.
    //  Cut the result according to (A.getNumCols())
    u.resize(t_A.getNumCols());
    std::string name;
    if (t_A.getName() != " " && t_v.getName() != " ") {
      name += "min_u||" + t_A.getName() + "*u - " + t_v.getName() + "|";
    }
    else {
      name  = " ";
    }
    Vector<double> U(name,u);
    if( info > 0 ) {
     std::string s = "The diagonal element " + std::to_string(info)
                   + " of the triangular factor of A is zero, so that A "
                   + "does not have full rank;the least squares solution"
                   + " could not be computed.";
     U.setFlag(-1);
     U.setInfo(s);
    }
    else if (info < 0) {
     U.setFlag(-1);
     U.setInfo(std::to_string(info)+"th argument has illegal value");
    }
    return U;
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DGELSS - solve the linear system min_C|B - AC|
  //           using SVD.
  //  Arguments:  A     - (m x k)-matrix
  //              B     - (m x n)-matrix
  //
  //  Returns:    Cm_min  ( (k x n)-matrix )
  //----------------------------------------------------------------------------
  Matrix<double> DGELSS(const Matrix<double>& t_A, const Matrix<double>& t_B)
  {
    if (t_A.getNumRows() != t_B.getNumRows()) {
      std::cout << "Matrices are incompatible!" << std::endl;
      return Matrix<double>("zeros",1,0.0);
    }
    Matrix<double> SVD(t_A);
    std::vector<double> c = t_B.getArray();
    int info;
    std::vector<double> singular(std::min(t_A.getNumRows(),t_A.getNumCols()));
    double rcond;
    int rank;
    info = LAPACKE_dgelss(LAPACK_ROW_MAJOR,  //  row major layout
                          t_A.getNumRows(),  //  number of rows of A
                          t_A.getNumCols(),  //  number of columns of A
                          t_B.getNumCols(),  //  number of columns of B
                          SVD.data(),        //  pointer to elements of A
                          t_A.getNumCols(),  //  leading dimension of A
                          c.data(),          //  pointer to elements of B
                          t_B.getNumCols(),  //  leading dimension of B
                          singular.data(),   //  pointer to the singular values
                          rcond,             //  condition number
                          &rank);            //  effective rank of A.
    //  Cut the result according to (A.getNumCols() x B.getNumCols())
    c.resize(t_A.getNumCols() * t_B.getNumCols());
    std::string name;
    if (t_A.getName() != " " && t_B.getName() != " ") {
      name += "min_C||" + t_A.getName() + "*C - " + t_B.getName() + "|";
    }
    else {
      name  = " ";
    }
    Matrix<double> C(name,t_A.getNumCols(),t_B.getNumCols(),c);
    if( info > 0 ) {
      std::string s = "The diagonal element " + std::to_string(info)
                    + " of the triangular factor of A is zero, so that A "
                    + "does not have full rank;the least squares solution"
                    + " could not be computed.";
      C.setFlag(-1);
      C.setInfo(s);
    }
    else if (info < 0) {
      C.setFlag(-1);
      C.setInfo(std::to_string(info)+"th argument has illegal value");
    }
    return C;
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DGELSS - solve the linear system min_u|v - Au|
  //           using SVD.
  //  Arguments:  A     - (m x k)-matrix
  //              v     - (m)-dim vector
  //
  //  Returns:    um_min  ( (k)-dim vector )
  //----------------------------------------------------------------------------
  Vector<double> DGELSS(const Matrix<double>& t_A, const Vector<double>& t_v)
  {
    if (t_A.getNumRows() != t_v.getDim()) {
      std::cout << "Matrices and vector are incompatible!" << std::endl;
      return Vector<double>("zeros",1,0.0);
    }
    Matrix<double> SVD(t_A);
    //  ------------------------------------------------------------------------
    //  WARNING: while A can be a general m x n matrix, the vector u
    //  must have at least as many rows as the number of columns of A,
    //  otherwise we get a segfault from LAPACK.  To prevent this, u
    //  is t_initialized with _dim = A.getNumCols() with zeros.  Then, the
    //  first v.getDim() elements are filled with v's elements.
    //  (N. Carrara - 6/30/2020)
    //  Error fixed with commit - bc555791c2ebc9aadebe44e5b74fbc367b7e2123.
    //--------------------------------------------------------------------------
    std::vector<double> u(std::max(t_A.getNumCols(),t_A.getNumRows()),0.0);
    for (auto i = 0; i < t_v.getDim(); i++)
    {
      u[i] = t_v(i);
    }
    int info;
    std::vector<double> singular(std::min(t_A.getNumRows(),t_A.getNumCols()));
    double rcond;
    int rank;
    info = LAPACKE_dgelss(LAPACK_ROW_MAJOR,  //  row major layout
                          t_A.getNumRows(),  //  number of rows of A
                          t_A.getNumCols(),  //  number of columns of A
                          1,                 //  dimension of v
                          SVD.data(),        //  pointer to elements of A
                          t_A.getNumCols(),  //  leading dimension of A
                          u.data(),          //  pointer to elements of v
                          1,                 //  leading dimension of v
                          singular.data(),   //  pointer to the singular values
                          rcond,             //  condition number
                          &rank);            //  effective rank of A.
    //  Cut the result according to (A.getNumCols())
    u.resize(t_A.getNumCols());
    std::string name;
    if (t_A.getName() != " " && t_v.getName() != " ") {
      name += "min_u||" + t_A.getName() + "*u - " + t_v.getName() + "|";
    }
    else {
      name  = " ";
    }
    Vector<double> U(name,u);
    if( info > 0 ) {
     std::string s = "The diagonal element " + std::to_string(info)
                   + " of the triangular factor of A is zero, so that A "
                   + "does not have full rank;the least squares solution"
                   + " could not be computed.";
     U.setFlag(-1);
     U.setInfo(s);
    }
    else if (info < 0) {
     U.setFlag(-1);
     U.setInfo(std::to_string(info)+"th argument has illegal value");
    }
    return U;
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DGETRF - construct an P^-1LU factorization
  //  Arguments:  A     - (m x m)-matrix
  //
  //  Returns:    std::vector<size_t> (list of pivots for P^-1)
  //----------------------------------------------------------------------------
  std::vector<size_t> DGETRF(const Matrix<double>& t_A)
  {
    if (t_A.getNumRows() != t_A.getNumCols()) {
      std::cout << "Matrix is not square!" << std::endl;
      return std::vector<size_t>(0);
    }
    Matrix<double> LU(t_A);
    std::vector<int> ipiv(std::min(t_A.getNumRows(),t_A.getNumCols()));
    int info;
    info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR,  // row major layout
                          t_A.getNumRows(),  // number of rows of A
                          t_A.getNumCols(),  // number of columns of A
                          LU.data(),         // pointer to the elements of A
                          t_A.getNumCols(),  // leading dimension of A
                          ipiv.data());      // pointer to pivot list
    std::vector<size_t> pivot(ipiv.begin(), ipiv.end());
    return pivot;
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DGETRF_L_U - compute an LU decomposition
  //  Arguments:  A     - (m x m)-matrix
  //
  //  Returns:    Matrix<double> (compound matrix of L and U)
  //              LU = [[u_11  u_12  ---  u_1m],
  //                    [l_21  u_22  ---  u_2m],
  //                    [ |     |     |    |  ],
  //                    [lm_m1  lm_m2  ---  um_mm]]
  //----------------------------------------------------------------------------
  Matrix<double> DGETRF_L_U(const Matrix<double>& t_A)
  {
    if (t_A.getNumRows() != t_A.getNumCols()) {
      std::cout << "Matrix is not square!" << std::endl;
      return Matrix<double>("zeros",1,0.0);
    }
    Matrix<double> LU(t_A);
    std::vector<int> ipiv(std::min(t_A.getNumRows(),t_A.getNumCols()));
    int info;
    info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR,  // row major layout
                          t_A.getNumRows(),  // number of rows of A
                          t_A.getNumCols(),  // number of columns of A
                          LU.data(),         // pointer to the elements of A
                          t_A.getNumCols(),  // leading dimension of A
                          ipiv.data());      // pointer to pivot list
    std::string name = "LU-Matrix";
    if (t_A.getName() != " ") {
      name += " of " + t_A.getName();
    }
    LU.setName(name);
    return LU;
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DGETRF_LU - compute an LU decomposition
  //  Arguments:  A     - (m x m)-matrix
  //
  //  Returns:    std::tuple<Matrix<double>,Matrix<double>> (L and U)
  //----------------------------------------------------------------------------
  std::tuple<Matrix<double>,Matrix<double>> DGETRF_LU(const Matrix<double>& t_A)
  {
    if (t_A.getNumRows() != t_A.getNumCols()) {
      std::cout << "Matrix is not square!" << std::endl;
      return {Matrix<double>("zeros",1,0.0),Matrix<double>("zeros",1,0.0)};
    }
    Matrix<double> LU(t_A);
    std::vector<int> ipiv(std::min(t_A.getNumRows(),t_A.getNumCols()));
    int info;
    info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR,  // row major layout
                          t_A.getNumRows(),  // number of rows of A
                          t_A.getNumCols(),  // number of columns of A
                          LU.data(),         // pointer to the elements of A
                          t_A.getNumCols(),  // leading dimension of A
                          ipiv.data());      // pointer to pivot list
    std::vector<double> l(t_A.getNumRows() * t_A.getNumRows(),0.0);
    std::vector<double> u(t_A.getNumRows() * t_A.getNumRows(),0.0);
    for (auto i = 0; i < t_A.getNumRows(); i++) {
      for (auto j = 0; j < t_A.getNumRows(); j++) {
        if (i == j) {
          l[i * t_A.getNumRows() + j] = 1.0;
          u[i * t_A.getNumRows() + j] = LU(i,j);
        }
        else if (j < i) {
          l[i * t_A.getNumRows() + j] = LU(i,j);
        }
        else {
          u[i * t_A.getNumRows() + j] = LU(i,j);
        }
      }
    }
    std::string name_l = "L-Matrix";
    std::string name_u = "U-Matrix";
    if (t_A.getName() != " ") {
      name_l += " of " + t_A.getName();
      name_u += " of " + t_A.getName();
    }
    Matrix<double> L(name_l,t_A.getNumRows(),l);
    Matrix<double> U(name_u,t_A.getNumRows(),u);
    return {L,U};
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DGETRF_PLU - compute an LU decomposition
  //  Arguments:  A     - (m x m)-matrix
  //
  //  Returns:    std::tuple<Matrix<double>,Matrix<double>,Matrix<double>>
  //              (P, L and U)
  //----------------------------------------------------------------------------
  std::tuple<Matrix<double>,Matrix<double>,Matrix<double>>
  DGETRF_PLU(const Matrix<double>& t_A)//  set the singular values of A
  {
    if (t_A.getNumRows() != t_A.getNumCols()) {
      std::cout << "Matrix is not square!" << std::endl;
      return {Matrix<double>("zeros",1,0.0),
              Matrix<double>("zeros",1,0.0),
              Matrix<double>("zeros",1,0.0)};
    }
    Matrix<double> LU(t_A);
    std::vector<int> ipiv(std::min(t_A.getNumRows(),t_A.getNumCols()));
    int info;
    info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR,  // row major layout
                          t_A.getNumRows(),  // number of rows of A
                          t_A.getNumCols(),  // number of columns of A
                          LU.data(),         // pointer to the elements of A
                          t_A.getNumCols(),  // leading dimension of A
                          ipiv.data());      // pointer to pivot list
    std::vector<double> l(t_A.getNumRows() * t_A.getNumRows(),0.0);
    std::vector<double> u(t_A.getNumRows() * t_A.getNumRows(),0.0);
    for (auto i = 0; i < t_A.getNumRows(); i++) {
      for (auto j = 0; j < t_A.getNumRows(); j++) {
        if (i == j) {
          l[i * t_A.getNumRows() + j] = 1.0;
          u[i * t_A.getNumRows() + j] = LU(i,j);
        }
        else if (j < i) {
          l[i * t_A.getNumRows() + j] = LU(i,j);
        }
        else {
          u[i * t_A.getNumRows() + j] = LU(i,j);
        }
      }
    }
    std::string name_p = "P-Matrix";
    std::string name_l = "L-Matrix";
    std::string name_u = "U-Matrix";
    if (t_A.getName() != " ") {
      name_p += " of " + t_A.getName();
      name_l += " of " + t_A.getName();
      name_u += " of " + t_A.getName();
    }
    std::vector<size_t> pivot(ipiv.begin(),ipiv.end());
    Matrix<double> P = permutationMatrix_d(t_A.getNumRows(),pivot);
    P.setName(name_p);
    Matrix<double> L(name_l,t_A.getNumRows(),l);
    Matrix<double> U(name_u,t_A.getNumRows(),u);
    return {P,L,U};
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DGEQRF - compute a QR decomposition
  //  Arguments:  A     - (m x n)-matrix
  //
  //  Returns:    Vector<double> (elementary reflectors)
  //----------------------------------------------------------------------------
  Vector<double> DGEQRF(const Matrix<double>& t_A)
  {
    Matrix<double> QR(t_A);
    std::vector<double> reflectors(std::min(t_A.getNumRows(),t_A.getNumCols()));
    int info;
    info = LAPACKE_dgeqrf(LAPACK_ROW_MAJOR,    // row major order
                          t_A.getNumRows(),    // number of rows of A
                          t_A.getNumCols(),    // number of columns of A
                          QR.data(),           // pointer to the elements of A
                          t_A.getNumCols(),    // leading dimension of A
                          reflectors.data());  // pointer to reflectors
    std::string name = "QR-reflectors";
    if (t_A.getName() != " ") {
      name += " of " + t_A.getName();
    }
    return Vector<double>(name,reflectors);
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DORGQF - compute a QR decomposition
  //  Arguments:  A     - (m x n)-matrix
  //              ref   - (min(m,n))-dim vector
  //
  //  Returns:    Matrix<double> (Q Matrix)
  //----------------------------------------------------------------------------
  Matrix<double> DORGQR(const Matrix<double>& t_A,
                             const Vector<double>& t_ref)
  {
    Matrix<double> QR(t_A);
    int info;
    Vector<double> reflect(t_ref);
    info = LAPACKE_dorgqr(LAPACK_ROW_MAJOR,  // row major order
                          t_A.getNumRows(),  // number of rows of A
                          t_A.getNumCols(),  // number of columns of A
                          t_ref.getDim(),    // number of elementary reflectors
                          QR.data(),         // pointer to the elements of A
                          t_A.getNumCols(),  // leading dimension of A
                          reflect.data());   // pointer to reflectors
    std::string name = "Q-Matrix";
    if (t_A.getName() != " ") {
      name += " of " + t_A.getName();
    }
    QR.setName(name);
    return QR;
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DGEQRF_QR - compute a QR decomposition
  //  Arguments:  A     - (m x n)-matrix
  //
  //  Returns:    std::tuple<Matrix<double>,Matrix<double>> (Q and R)
  //----------------------------------------------------------------------------
  std::tuple<Matrix<double>,Matrix<double>>
  DGEQRF_QR(const Matrix<double>& t_A)
  {
    Matrix<double> Q(t_A);
    std::vector<double> reflectors(std::min(t_A.getNumRows(),t_A.getNumCols()));
    int info;
    info = LAPACKE_dgeqrf(LAPACK_ROW_MAJOR,    // row major order
                          t_A.getNumRows(),    // number of rows of A
                          t_A.getNumCols(),    // number of columns of A
                          Q.data(),            // pointer to the elements of A
                          t_A.getNumCols(),    // leading dimension of A
                          reflectors.data());  // pointer to reflectors
    //  find Q using dorgqr
    info = LAPACKE_dorgqr(LAPACK_ROW_MAJOR,    // row major order
                          t_A.getNumRows(),    // number of rows of A
                          t_A.getNumCols(),    // number of columns of A
                          reflectors.size(),   // number of elementary ref
                          Q.data(),            // pointer to the elements of A
                          t_A.getNumCols(),    // leading dimension of A
                          reflectors.data());  // pointer to reflectors
    std::string name_q = "Q-Matrix";
    std::string name_r = "R-Matrix";
    if (t_A.getName() != " ") {
      name_q += " of " + t_A.getName();
      name_r += " of " + t_A.getName();
    }
    Q.setName(name_q);
    Matrix<double> Q_T(Q);
    Q_T.transpose_inplace();
    //  find R = Q^T*A
    Matrix<double> R = Q_T * t_A;
    R.setName(name_r);
    return {Q,R};
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DGESVD - compute a SVD decomposition
  //  Arguments:  A     - (m x n)-matrix
  //
  //  Returns:    Vector<double> (singular values of A)
  //----------------------------------------------------------------------------
  Vector<double> DGESVD(const Matrix<double>& t_A)
  {
    Matrix<double> A_copy(t_A);
    Matrix<double> U(t_A.getNumRows());
    Matrix<double> VT(t_A.getNumCols());
    double superb[std::min(t_A.getNumRows(),t_A.getNumCols())-1];
    std::vector<double> singular(std::min(t_A.getNumRows(),t_A.getNumCols()));
    int info;
    info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR,// row major format
                          'N',             // no columns returned in U
                          'N',             // no rows returned in VT
                          t_A.getNumRows(),// number of rows of A
                          t_A.getNumCols(),// number of columns of A
                          A_copy.data(),   // pointer to the elements of A
                          t_A.getNumCols(),// leading dimension of A
                          singular.data(), // pointer to singular values
                          U.data(),        // pointer to elements of U
                          U.getNumCols(),  // leading dimension of U
                          VT.data(),       // pointer to elements of VT
                          VT.getNumCols(), // leading dimension of VT
                          superb);         // workspace for subroutines
    std::string name = "Singular values";
    if (t_A.getName() != " ") {
      name += " of " + t_A.getName();
    }
    return Vector<double>(name,singular);
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DGESVD - compute a SVD decomposition
  //  Arguments:  A     - (m x n)-matrix
  //
  //  Returns:    std::tuple<Matrix<double>,Matrix<double>,Matrix<double>>
  //              (Matrices U, Sigma and (V)^T)
  //----------------------------------------------------------------------------
  std::tuple<Matrix<double>,Matrix<double>,Matrix<double>>
  DGESVD_SVD(const Matrix<double>& t_A)
  {
    Matrix<double> A_copy(t_A);
    Matrix<double> U(t_A.getNumRows());
    Matrix<double> VT(t_A.getNumCols());
    double superb[std::min(t_A.getNumRows(),t_A.getNumCols())-1];
    std::vector<double> singular(std::min(t_A.getNumRows(),t_A.getNumCols()));
    int info;
    info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR,// row major format
                          'A',             // all m rows returned in U
                          'A',             // all n rows returned in VT
                          t_A.getNumRows(),// number of rows of A
                          t_A.getNumCols(),// number of columns of A
                          A_copy.data(),   // pointer to the elements of A
                          t_A.getNumCols(),// leading dimension of A
                          singular.data(), // pointer to singular values
                          U.data(),        // pointer to elements of U
                          U.getNumCols(),  // leading dimension of U
                          VT.data(),       // pointer to elements of VT
                          VT.getNumCols(), // leading dimension of VT
                          superb);         // workspace for subroutines
    std::string name_u = "U-Matrix";
    std::string name_s = "Sigma-Matrix";
    std::string name_vt = "(V)^T-Matrix";
    if (t_A.getName() != " ") {
      name_u += " of " + t_A.getName();
      name_s += " of " + t_A.getName();
      name_vt += " of " + t_A.getName();
    }
    U.setName(name_u);
    VT.setName(name_vt);
    Matrix<double> Sigma(name_s,t_A.getNumRows(),t_A.getNumCols(),0.0);
    for (auto i = 0; i < singular.size(); i++) {
      Sigma(i,i) = singular[i];
    }
    return {U,Sigma,VT};
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DGESDD - compute a SVD decomposition
  //  Arguments:  A     - (m x n)-matrix
  //
  //  Returns:    Vector<double> (singular values of A)
  //----------------------------------------------------------------------------
  Vector<double> DGESDD(const Matrix<double>& t_A)
  {
    Matrix<double> A_copy(t_A);
    Matrix<double> U(t_A.getNumRows());
    Matrix<double> VT(t_A.getNumCols());
    double superb[std::min(t_A.getNumRows(),t_A.getNumCols())-1];
    std::vector<double> singular(std::min(t_A.getNumRows(),t_A.getNumCols()));
    int info;
    info = LAPACKE_dgesdd(LAPACK_ROW_MAJOR,// row major format
                          'N',             // no columns returned in U or V
                          t_A.getNumRows(),// number of rows of A
                          t_A.getNumCols(),// number of columns of A
                          A_copy.data(),   // pointer to the elements of A
                          t_A.getNumCols(),// leading dimension of A
                          singular.data(), // pointer to singular values
                          U.data(),        // pointer to elements of U
                          U.getNumCols(),  // leading dimension of U
                          VT.data(),       // pointer to elements of VT
                          VT.getNumCols());// leading dimension of VT
    std::string name = "Singular values";
    if (t_A.getName() != " ") {
      name += " of " + t_A.getName();
    }
    return Vector<double>(name,singular);
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DGESDD - compute a SVD decomposition
  //  Arguments:  A     - (m x n)-matrix
  //
  //  Returns:    std::tuple<Matrix<double>,Matrix<double>,Matrix<double>>
  //              (Matrices U, Sigma and (V)^T)
  //----------------------------------------------------------------------------
  std::tuple<Matrix<double>,Matrix<double>,Matrix<double>>
  DGESDD_SVD(const Matrix<double>& t_A)
  {
    Matrix<double> A_copy(t_A);
    Matrix<double> U(t_A.getNumRows());
    Matrix<double> VT(t_A.getNumCols());
    double superb[std::min(t_A.getNumRows(),t_A.getNumCols())-1];
    std::vector<double> singular(std::min(t_A.getNumRows(),t_A.getNumCols()));
    int info;
    info = LAPACKE_dgesdd(LAPACK_ROW_MAJOR,// row major format
                          'A',             // all m rows returned in U and VT
                          t_A.getNumRows(),// number of rows of A
                          t_A.getNumCols(),// number of columns of A
                          A_copy.data(),   // pointer to the elements of A
                          t_A.getNumCols(),// leading dimension of A
                          singular.data(), // pointer to singular values
                          U.data(),        // pointer to elements of U
                          U.getNumCols(),  // leading dimension of U
                          VT.data(),       // pointer to elements of VT
                          VT.getNumCols());// leading dimension of VT
    std::string name_u = "U-Matrix";
    std::string name_s = "Sigma-Matrix";
    std::string name_vt = "(V)^T-Matrix";
    if (t_A.getName() != " ") {
      name_u += " of " + t_A.getName();
      name_s += " of " + t_A.getName();
      name_vt += " of " + t_A.getName();
    }
    U.setName(name_u);
    VT.setName(name_vt);
    Matrix<double> Sigma(name_s,t_A.getNumRows(),t_A.getNumCols(),0.0);
    for (auto i = 0; i < singular.size(); i++) {
      Sigma(i,i) = singular[i];
    }
    return {U,Sigma,VT};
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DGETRI - compute the inverse of a matrix using the LU decomposition
  //  Arguments:  A     - (m x m)-matrix
  //
  //  Returns:    Matrix<T> inverse
  //----------------------------------------------------------------------------
  Matrix<double> DGETRI(const Matrix<double>& t_A)
  {
    if (t_A.getNumRows() != t_A.getNumCols()) {
      std::cout << "Matrix is not square!" << std::endl;
      return t_A;
    }
    Matrix<double> A_copy(t_A);
    Matrix<double> LU(t_A);
    //  First compute the LU factorization to get the
    //  pivot indices
    std::vector<int> ipiv(std::min(t_A.getNumRows(),t_A.getNumCols()));
    int info;
    info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR,// row major layout
                          t_A.getNumRows(),// number of rows of A
                          t_A.getNumCols(),// number of columns of A
                          LU.data(),       // pointer to the elements of A
                          t_A.getNumCols(),// leading dimension of A
                          ipiv.data());    // pointer to pivot list
    //  Now to take the inverse
    info = LAPACKE_dgetri(LAPACK_ROW_MAJOR,// row major layout
                          t_A.getNumRows(),// number of rows of A
                          LU.data(),       // pointer to the elements of A
                          t_A.getNumCols(),// leading dimension of A
                          ipiv.data());    // pointer to the pivot list
    if (info > 0) {
      LU.setFlag(-1);
      LU.setInfo("Matrix could not be computed, diagonal element "
                 + std::to_string(info) + " is exactly zero");
    }
    std::string name;
    if (t_A.getName() != " ") {
      name += "(" + t_A.getName() + ")^-1";
    }
    LU.setName(name);
    return LU;
  }
  //----------------------------------------------------------------------------
}
