//------------------------------------------------------------------------------
//  matrix.cpp
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
#include "vector.h"
#include "matrix.h"

//------------------------------------------------------------------------------
//  matrix.cpp
//
//  This file contains the instantiations of the definitions given in matrix.h.
//  The various methods are listed below,
//  ----------------------------------------------------------------------------
//  Line no.  |   Method
//  ----------------------------------------------------------------------------
//  (114)     |   Matrix()
//  (121)     |   ~Matrix()
//  (129)     |   Matrix(const Matrix<T>&)
//  (149)     |   Matrix(uint32_t)
//  (154)     |   Matrix(std::string, uint32_t)
//  (160)     |   Matrix(uint32_t, uint32_t)
//  (165)     |   Matrix(std::string, uint32_t, uint32_t)
//  (171)     |   Matrix(uint32_t, uint32_t, const T&)
//  (177)     |   Matrix(std::string, uint32_t, uint32_t, const T&)
//  (186)     |   Matrix(uint32_t, std::vector<T>)
//  (191)     |   Matrix(std::string, uint32_t, std::vector<T>)
//  (196)     |   Matrix(uint32_t, uint32_t, std::vector<T>)
//  (201)     |   Matrix(std::string, uint32_t, uint32_t, std::vector<T>)
//  (206)     |   Matrix(std::string, uint32_t, uint32_t, T*)
//  (213)     |   Matrix(std::vector<std::vector<T>>)
//  (224)     |   Matrix(std::string, std::vector<std::vector<T>>)
//  (243)     |   getNumRows()
//  (248)     |   getNumCols()
//  (253)     |   getName()
//  (259)     |   getArray()
//  (267)     |   accessArray()
//  (273)     |   data()
//  (278)     |   getRow(uint32_t)
//  (295)     |   getCol(uint32_t)
//  (306)     |   setName(std::string)
//  (311)     |   setRow(uint32_t, std::vector<T>)
//  (319)     |   setCol(uint32_t, std::vector<T>)
//  (327)     |   setArray(uint32_t, std::vector<T>)
//  (334)     |   setArray(std::vector<std::vector<T>>)
//  (353)     |   operator=(const Matrix<T>&)
//  (370)     |   operator==(const Matrix<T>&)
//  (386)     |   operator!=(const Matrix<T>&)
//  (402)     |   operator-()
//  (413)     |   operator+(const Matrix<T>&)
//  (429)     |   operator+=(const Matrix<T>&)
//  (445)     |   operator-(const Matrix<T>&)
//  (461)     |   operator-=(const Matrix<T>&)
//  (477)     |   operator*(const Matrix<T>&)
//  (488)     |   brute_mul(const Matrix<T>&)
//  (507)     |   operator*=(const Matrix<T>&)
//  (515)     |   operator+(const T&)
//  (526)     |   operator-(const T&)
//  (537)     |   operator*(const T&)
//  (548)     |   operator/(const T&)
//  (564)     |   operator+=(const T&)
//  (575)     |   operator-=(const T&)
//  (586)     |   operator*=(const T&)
//  (597)     |   operator/=(const T&)
//  (613)     |   operator*(const Vector<T>&)
//  (634)     |   operator()(const uint32_t&, const uint32_t&)
//  (639)     |   operator()(const uint32_t&, const uint32_t&)
//  (644)     |   operator()(const uint32_t&)
//  (649)     |   operator()(const uint32_t&)
//  (664)     |   print()
//  (682)     |   summary()
//  (882)     |   transpose()
//  (898)     |   transpose_inplace(bool)
//  (917)     |   trace()
//  ()        |   permutationMatrix()
//  ()        |   getRank()
//  ()        |   isInvertible()
//  ()        |   findSingularValues()
//  ()        |   inverse()
//  ()        |   pseudoInverse()
//  ()        |   LU()
//  ()        |   getL(const Matrix<T>&)
//  ()        |   getU(const Matrix<T>&)
//  ()        |   QR()
//  ()        |   SVD()
//  ()        |   getSingularValues()
//  ()        |   identity_d(uint32_t)
//------------------------------------------------------------------------------

namespace ET
{
  //----------------------------------------------------------------------------
  //  Matrix constructors
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Default constructor
  //    sets name = " ", and _m, _n = 0
  //----------------------------------------------------------------------------
  template<typename T>
  Matrix<T>::Matrix() : _name(" "), _m(0), _n(0)
  {
    std::cout << "\nMatrix created at location " << this;
  }
  //----------------------------------------------------------------------------
  //  Default destructor
  //----------------------------------------------------------------------------
  template<typename T>
  Matrix<T>::~Matrix()
  {
    //std::cout << "\nMatrix " + _name + " at location " << this << " destroyed.";
  }
  //----------------------------------------------------------------------------
  //  Copy constructor
  //    Does not delete the copied object.
  //----------------------------------------------------------------------------
  template<typename T>
  Matrix<T>::Matrix(const Matrix<T>& matrix)
  {
    //std::cout << "\nMatrix created at location " << this;
    _array = matrix._array;
    _m = matrix.getNumRows();
    _n = matrix.getNumCols();
    _name = matrix.getName();
  }
  //----------------------------------------------------------------------------
  //  Constructors with various sets of arguments, such as,
  //    std::string                 name,
  //    uint32_t                    m,
  //    uint32_t                    n,
  //    std::vector<T>              flatted array,
  //    std::vector<std::vector<T>> 2d array,
  //    T*                          C style flattened array,
  //    const T&                    initial value for all elements.
  //  If any sizes are specified, the internal array element '_array'
  //  is resized accordingly, otherwise it is left uninstantiated.
  //----------------------------------------------------------------------------
  template<typename T>
  Matrix<T>::Matrix(uint32_t m) : _m(m), _n(m), _name(" ")
  {
    //std::cout << "\nMatrix created at location " << this;
    _array.resize(_m*_n,0.0);
  }
  template<typename T>
  Matrix<T>::Matrix(std::string name, uint32_t m)
  : _m(m), _n(m), _name(name)
  {
    //std::cout << "\nMatrix created at location " << this;
    _array.resize(_m*_n,0.0);
  }
  template<typename T>
  Matrix<T>::Matrix(uint32_t m, uint32_t n) : _m(m), _n(n), _name(" ")
  {
    //std::cout << "\nMatrix created at location " << this;
    _array.resize(_m*_n,0.0);
  }
  template<typename T>
  Matrix<T>::Matrix(std::string name, uint32_t m, uint32_t n)
  : _m(m), _n(n), _name(name)
  {
    //std::cout << "\nMatrix created at location " << this;
    _array.resize(_m*_n,0.0);
  }
  template<typename T>
  Matrix<T>::Matrix(uint32_t m, uint32_t n, const T& init)
  : _m(m), _n(n), _name(" ")
  {
    //std::cout << "\nMatrix created at location " << this;
    _array.resize(_m*_n, init);
  }
  template<typename T>
  Matrix<T>::Matrix(std::string name, uint32_t m,
    uint32_t n, const T& init)
  : _m(m), _n(n), _name(name)
  {
    //std::cout << "\nMatrix created at location " << this;
    _array.resize(_m*_n, init);
  }
  //  Notice that the following methods do not MOVE the vectors so that they
  //  change ownership.  Instead they are copied into _array.
  template<typename T>
  Matrix<T>::Matrix(uint32_t m, std::vector<T> flat)
  : _m(m), _n(m), _name(" "), _array(flat)
  {
    std::cout << "\nMatrix created at location " << this;
  }
  template<typename T>
  Matrix<T>::Matrix(std::string name, uint32_t m, std::vector<T> flat)
  : _m(m), _n(m), _name(name), _array(flat)
  {
    //std::cout << "\nMatrix created at location " << this;
  }
  template<typename T>
  Matrix<T>::Matrix(uint32_t m, uint32_t n, std::vector<T> flat)
  : _m(m), _n(n), _name(" "), _array(flat)
  {
    //std::cout << "\nMatrix created at location " << this;
  }
  template<typename T>
  Matrix<T>::Matrix(std::string name, uint32_t m, uint32_t n, std::vector<T> flat)
  : _m(m), _n(n), _name(name), _array(flat)
  {
    //std::cout << "\nMatrix created at location " << this;
  }
  template<typename T>
  Matrix<T>::Matrix(std::string name, uint32_t m, uint32_t n, T* array)
  : _m(m), _n(n), _name(name)
  {
    //std::cout << "\nMatrix created at location " << this;
    std::vector<T> flat(array, array + _m*_n);
    _array = flat;
  }
  template<typename T>
  Matrix<T>::Matrix(std::vector<std::vector<T> > array)
  : _m(array.size()), _n(array[0].size()), _name(" ")
  {
    //std::cout << "\nMatrix created at location " << this;
    std::vector<T> flat;
    for (uint32_t i = 0; i < _m; i++)
    {
      flat.insert(end(flat),begin(array[i]),end(array[i]));
    }
    _array = flat;
  }
  template<typename T>
  Matrix<T>::Matrix(std::string name, std::vector<std::vector<T> > array)
  : _m(array.size()), _n(array[0].size()), _name(name)
  {
    //std::cout << "\nMatrix created at location " << this;
    std::vector<std::pair<uint32_t,uint32_t>> rows;
    bool well_defined = checkConsistency(array,rows);
    if (well_defined == false)
    {
      setInfo(MATRIX_INCONSISTENT_ARRAY(rows));
    }
    std::vector<T> flat;
    for (uint32_t i = 0; i < _m; i++)
    {
      flat.insert(end(flat),begin(array[i]),end(array[i]));
    }
    _array = flat;
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Getters and Setters
  //  Each attribute comes with its own setters and getters.  There are a few
  //  additional ones such as 'getRow' and 'getCol' which each return an
  //  std::vector<T>.
  //----------------------------------------------------------------------------
  template<typename T>
  uint32_t Matrix<T>::getNumRows() const
  {
    return _m;
  }
  template<typename T>
  uint32_t Matrix<T>::getNumCols() const
  {
    return _n;
  }
  template<typename T>
  std::string Matrix<T>::getName() const
  {
    return _name;
  }
  //  When 'getArray()' is called it will usually create a copy of _array.
  template<typename T>
  std::vector<T> Matrix<T>::getArray() const
  {
    return _array;
  }
  //  In order to return the '_array' attribute so that it can be manipulated
  //  by other methods, such as those which utilize BLAS and LAPACK functions,
  //  we use 'accessArray()' to return a pointer to '_array'.
  template<typename T>
  std::vector<T>* Matrix<T>::accessArray()
  {
    return &_array;
  }
  //  get access to the beginning of the array
  template<typename T>
  T* Matrix<T>::data()
  {
    return _array.data();
  }
  template<typename T>
  std::vector<T> Matrix<T>::getRow(uint32_t i)
  {
    std::vector<T> row(_m,0.0);
    //  check that row exists, if not log error and send zero vector
    if (i >= _m)
    {
      _info = MATRIX_OUT_OF_BOUNDS(0,_n,i,_name);
      return row;
    }
    for (uint32_t j = 0; j < _m; j++)
    {
      row[j] = _array[i*_m + j];
    }
    return row;
  }
  template<typename T>
  std::vector<T> Matrix<T>::getCol(uint32_t i)
  {
    std::vector<T> col(_m,0.0);
    for (uint32_t j = 0; j < _n; j++)
    {
      col[j] = _array[j*_n + i];
    }
    return col;
  }
  template<typename T>
  std::vector<T> Matrix<T>::getSingularValues()
  {
    return _singular_values;
  }
  template<typename T>
  std::string Matrix<T>::getInfo()
  {
    return _info;
  }
  template<typename T>
  int Matrix<T>::getFlag()
  {
    return _flag;
  }
  template<typename T>
  uint32_t Matrix<T>::getRank()
  {
    return _rank;
  }
  //  Setters
  template<typename T>
  void Matrix<T>::setName(std::string name)
  {
    _name = name;
  }
  template<typename T>
  void Matrix<T>::setRow(uint32_t i, std::vector<T> row)
  {
    for (uint32_t j = 0; j < _n; j++)
    {
      _array[i*_n + j] = row[j];
    }
  }
  template<typename T>
  void Matrix<T>::setCol(uint32_t i, std::vector<T> col)
  {
    for (uint32_t j = 0; j < _m; j++)
    {
      _array[j*_n + i] = col[j];
    }
  }
  template<typename T>
  void Matrix<T>::setArray(uint32_t m, std::vector<T> mat)
  {
    _n = mat.size()/m;
    _m = m;
    _array = mat;
  }
  template<typename T>
  void Matrix<T>::setArray(std::vector<std::vector<T> > mat)
  {
    std::vector<std::pair<uint32_t,uint32_t>> rows;
    bool well_defined = checkConsistency(mat,rows);
    if (well_defined == false)
    {
      setInfo(MATRIX_INCONSISTENT_ARRAY(rows));
    }
    _m = mat.size();
    _n = mat[0].size();
    _array.resize(_m*_n);
    for (uint32_t i = 0; i < _m; i++)
    {
      for (uint32_t j = 0; j < _n; j++)
      {
        _array[i*_n + j] = mat[i][j];
      }
    }
  }
  template<typename T>
  void Matrix<T>::setSingularValues(std::vector<T> singular)
  {
    _singular_values = singular;
  }
  template<typename T>
  void Matrix<T>::setInfo(std::string info)
  {
    _info = info;
  }
  template<typename T>
  void Matrix<T>::setFlag(int flag)
  {
    _flag = flag;
  }
  template<typename T>
  void Matrix<T>::setRank(uint32_t rank)
  {
    _rank = rank;
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Operator overloads
  //----------------------------------------------------------------------------
  template<typename T>
  Matrix<T>& Matrix<T>::operator=(const Matrix<T>& matrix)
  {
    if (&matrix == this)
    {
      return *this;
    }
    _m = matrix.getNumRows();
    _n = matrix.getNumCols();
    _name = matrix.getName();
    _array.resize(_m*_n);
    for (uint32_t i = 0; i < _m*_n; i++)
    {
        _array[i] = matrix(i);
    }
    return *this;
  }
  template<typename T>
  bool Matrix<T>::operator==(const Matrix<T>& matrix) const
  {
    if (_n != matrix.getNumRows() || _m != matrix.getNumCols())
    {
      return false;
    }
    for (uint32_t i = 0; i < _m*_n; i++)
    {
        if (matrix(i) != _array[i])
        {
          return false;
        }
    }
    return true;
  }
  template<typename T>
  bool Matrix<T>::operator!=(const Matrix<T>& matrix) const
  {
    if (_n != matrix.getNumRows() || _m != matrix.getNumCols())
    {
      return true;
    }
    for (uint32_t i = 0; i < _m*_n; i++)
    {
        if (matrix(i) != _array[i])
        {
          return true;
        }
    }
    return false;
  }
  template<typename T>
  Matrix<T> Matrix<T>::operator-() const
  {
    std::vector<T> mat(_m*_n);
    for (uint32_t i = 0; i < _m*_n; i++)
    {
      mat[i] = -1*_array[i];
    }
    return Matrix<T>(_name,_m,_n,mat);
  }
  //  Matrix algebra
  template<typename T>
  Matrix<T> Matrix<T>::operator+(const Matrix<T>& matrix) const
  {
    if(_n != matrix.getNumRows())
    {
      setInfo(MATRIX_ADD_INCOMPATIBLE_ROWS(_n, matrix.getNumRows(),
                                           _name, matrix.getName()));
      return *this;
    }
    if (_m != matrix.getNumCols())
    {
      setInfo(MATRIX_ADD_INCOMPATIBLE_COLS(_m, matrix.getNumCols(),
                                           _name, matrix.getName()));
      return *this;
    }
    std::string name = "(" + _name + " + " + matrix.getName() + ")";
    Matrix<T> l(name, _m, _n, 0.0);
    for (uint32_t i = 0; i < _m*_n; i++)
    {
        l(i) = _array[i] + matrix(i);
    }
    return l;
  }
  template<typename T>
  Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& matrix)
  {
    if(_n != matrix.getNumRows() || _m != matrix.getNumCols())
    {
      std::cout << "Matrices incompatible!" << std::endl;
      return *this;
    }
    std::string name = "(" + _name + " + " + matrix.getName() + ")";
    setName(name);
    for (uint32_t i = 0; i < _m*_n; i++)
    {
        _array[i] += matrix(i);
    }
    return *this;
  }
  template<typename T>
  Matrix<T> Matrix<T>::operator-(const Matrix<T>& matrix) const
  {
    if(_n != matrix.getNumRows() || _m != matrix.getNumCols())
    {
      std::cout << "Matrices incompatible!" << std::endl;
      return *this;
    }
    std::string name = "(" + _name + " - " + matrix.getName() + ")";
    Matrix<T> l(name,_m,_n,0.0);
    for (uint32_t i = 0; i < _m*_n; i++)
    {
        l(i) = _array[i] - matrix(i);
    }
    return l;
  }
  template<typename T>
  Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& matrix)
  {
    if(_n != matrix.getNumRows() || _m != matrix.getNumCols())
    {
      std::cout << "Matrices incompatible!" << std::endl;
      return *this;
    }
    std::string name = "(" + _name + " - " + matrix.getName() + ")";
    setName(name);
    for (uint32_t i = 0; i < _m*_n; i++)
    {
        _array[i] -= matrix(i);
    }
    return *this;
  }
  template<typename T>
  Matrix<T> Matrix<T>::operator*(const Matrix<T>& matrix) const
  {
    if(_n != matrix.getNumRows())
    {
      std::cout << "Matrices " + _name + " and " + matrix.getName()
          + " incompatible!" << std::endl;
      return *this;
    }
    return DGEMM(*this,matrix);
  }
  template<typename T>
  Matrix<T> Matrix<T>::brute_mul(const Matrix<T>& matrix) const
  {
    if(_n != matrix.getNumRows())
    {
      std::cout << "Matrices incompatible!" << std::endl;
      return *this;
    }
    std::string name = "(" + _name + " * " + matrix.getName() + ")";
    Matrix<T> l(name,_m,matrix.getNumCols(),0.0);
    for (uint32_t i = 0; i < _m; i++) {
      for (uint32_t j = 0; j < _n; j++) {
        for (uint32_t k = 0; k < _m; k++) {
          l(i,j) += this->_array[i*_n + k] * matrix(k,j);
        }
      }
    }
    return l;
  }
  template<typename T>
  Matrix<T>& Matrix<T>::operator*=(const Matrix<T>& matrix)
  {
    Matrix<T> l = (*this) * matrix;
    (*this) = l;
    return *this;
  }
  //  Scalar algebra
  template<typename T>
  Matrix<T> Matrix<T>::operator+(const T& s) const
  {
    std::string name = "(" + _name + " + " + std::to_string(s) + "I)";
    Matrix<T> l(name,_m,_n,0.0);
    for (uint32_t i = 0; i < _m*_n; i++)
    {
        l(i) = _array[i] + s;
    }
    return l;
  }
  template<typename T>
  Matrix<T> Matrix<T>::operator-(const T& s) const
  {
    std::string name = "(" + _name + " - " + std::to_string(s) + "I)";
    Matrix<T> l(name,_m,_n,0.0);
    for (uint32_t i = 0; i < _m*_n; i++)
    {
        l(i) = _array[i] - s;
    }
    return l;
  }
  template<typename T>
  Matrix<T> Matrix<T>::operator*(const T& s) const
  {
    std::string name = "(" + _name + " * " + std::to_string(s) + ")";
    Matrix<T> l(name,_m,_n,0.0);
    for (uint32_t i = 0; i < _m*_n; i++)
    {
        l(i) = _array[i] * s;
    }
    return l;
  }
  template<typename T>
  Matrix<T> Matrix<T>::operator/(const T& s) const
  {
    std::string name = "(" + _name + " / " + std::to_string(s) + ")";
    Matrix<T> l(name,_m,_n,0.0);
    if (s == 0)
    {
      return l;
    }
    for (uint32_t i = 0; i < _m*_n; i++)
    {
        l(i) = _array[i] + s;
    }
    return l;
  }

  template<typename T>
  Matrix<T>& Matrix<T>::operator+=(const T& s)
  {
    std::string name = "(" + _name + " + " + std::to_string(s) + "I)";
    setName(name);
    for (uint32_t i = 0; i < _m*_n; i++)
    {
        _array[i] += s;
    }
    return *this;
  }
  template<typename T>
  Matrix<T>& Matrix<T>::operator-=(const T& s)
  {
    std::string name = "(" + _name + " - " + std::to_string(s) + "I)";
    setName(name);
    for (uint32_t i = 0; i < _m*_n; i++)
    {
        _array[i] -= s;
    }
    return *this;
  }
  template<typename T>
  Matrix<T>& Matrix<T>::operator*=(const T& s)
  {
    std::string name = "(" + _name + " * " + std::to_string(s) + ")";
    setName(name);
    for (uint32_t i = 0; i < _m*_n; i++)
    {
        _array[i] *= s;
    }
    return *this;
  }
  template<typename T>
  Matrix<T>& Matrix<T>::operator/=(const T& s)
  {
    if (s == 0)
    {
      return *this;
    }
    std::string name = "(" + _name + " / " + std::to_string(s) + ")";
    setName(name);
    for (uint32_t i = 0; i < _m*_n; i++)
    {
        _array[i] += s;
    }
    return *this;
  }
  //  Multiplying a vector
  template<typename T>
  Vector<T> Matrix<T>::operator*(const Vector<T>& v)
  {
    std::vector<T> vec(_n,0.0);
    Vector<T> v2(vec);
    if (_n != v.getDim())
    {
      return v2;
    }
    for (uint32_t i = 0; i < _m; i++)
    {
      T temp = 0.0;
      for (uint32_t j = 0; j < _n; j++)
      {
        temp += this->_array[i*_n + j] * v(j);
      }
      v2(i) = temp;
    }
    return v2;
  }
  //  Access Matrix<T>::operators
  template<typename T>
  T& Matrix<T>::operator()(const uint32_t& i, const uint32_t& j)
  {
    return _array[i*_n + j];
  }
  template<typename T>
  const T& Matrix<T>::operator()(const uint32_t& i, const uint32_t& j) const
  {
    return _array[i*_n + j];
  }
  template<typename T>
  T& Matrix<T>::operator()(const uint32_t& i)
  {
    return _array[i];
  }
  template<typename T>
  const T& Matrix<T>::operator()(const uint32_t& i) const
  {
    return _array[i];
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
    std::cout << "(" << _m << " x " << _n << ") Matrix";
    if (_name != " ")
      std::cout << ": '" << _name << "'";
    std:: cout << "\n[ ";
    for (uint32_t i = 0; i < _n; i++) {
      for (uint32_t j = 0; j < _m; j++) {
        std::cout << this->_array[i*_m + j];
        if (j < _m-1)
        std::cout << "   ";
      }
      if (i < _n-1)
        std::cout << "\n  ";
    }
    std::cout << " ]" << std::endl;
  }
  template<typename T>
  std::string Matrix<T>::summary()
  {
    std::stringstream s;
    s.str("");
    s.clear();
    std::string sum = "dim: (" + std::to_string(_m)
                    + "x" + std::to_string(_n)
                    + "), type: "
                    + type_name<decltype(_array[0])>();
    if (_name != " ")
    {
      sum +=  ", name: '" + _name + "'";
    }
    if (_array.size() == 0)
    {
      sum += "\n[  empty  ]";
      return sum;
    }
    sum += "\n[ ";
    if (_m < 10)
    {
      for (uint32_t i = 0; i < _m; i++)
      {
        if (_n < 10)
        {
          if (_array[i*_n] >= 0.0)
            sum += " ";
          for (uint32_t j = 0; j < _n; j++)
          {
            sum += scientific_not(this->_array[i*_n + j],3);
            if (j < _n-1)
            {
              if (_array[i*_n + j+1] >= 0.0)
                sum += "   ";
              else
                sum += "  ";
            }
          }
        }
        else
        {
          if (_array[i*_n] >= 0.0)
            sum += " ";
          sum += scientific_not(this->_array[i*_n + 0],3);
          if (_array[i*_n + 1] >= 0.0)
            sum += "   ";
          else
            sum += "  ";
          sum += scientific_not(this->_array[i*_n + 1],3);
          if (_array[i*_n + 2] >= 0.0)
            sum += "   ";
          else
            sum += "  ";
          sum += scientific_not(this->_array[i*_n + 2],3);
          sum += "   ";
          sum += "...   ...   ...   ";
          if (_array[i*_n + _n-3] >= 0.0)
            sum += " ";
          sum += scientific_not(this->_array[i*_n + _n-3],3);
          if (_array[i*_n + _n-2] >= 0.0)
            sum += "   ";
          else
            sum += "  ";
          sum += scientific_not(this->_array[i*_n + _n-2],3);
          if (_array[i*_n + _n-1] >= 0.0)
            sum += "   ";
          else
            sum += "  ";
          sum += scientific_not(this->_array[i*_n + _n-1],3);
        }
        if (i < _m-1)
        {
          sum += "\n  ";
        }
      }
    }
    else
    {
      for (uint32_t i = 0; i < 3; i++)
      {
        if (_n < 10)
        {
          if (_array[i*_n] >= 0.0)
            sum += " ";
          for (uint32_t j = 0; j < _n; j++)
          {
            sum += scientific_not(this->_array[i*_n + j],3);
            if (j < _n-1)
            {
              if (_array[i*_n + j+1] >= 0.0)
                sum += "   ";
              else
                sum += "  ";
            }
          }
        }
        else
        {
          if (_array[i*_n] >= 0.0)
            sum += " ";
          sum += scientific_not(this->_array[i*_n + 0],3);
          if (_array[i*_n + 1] >= 0.0)
            sum += "   ";
          else
            sum += "  ";
          sum += scientific_not(this->_array[i*_n + 1],3);
          if (_array[i*_n + 2] >= 0.0)
            sum += "   ";
          else
            sum += "  ";
          sum += scientific_not(this->_array[i*_n + 2],3);
          sum += "   ";
          sum += "...   ...   ...   ";
          if (_array[i*_n + _n-3] >= 0.0)
            sum += " ";
          sum += scientific_not(this->_array[i*_n + _n-3],3);
          if (_array[i*_n + _n-2] >= 0.0)
            sum += "   ";
          else
            sum += "  ";
          sum += scientific_not(this->_array[i*_n + _n-2],3);
          if (_array[i*_n + _n-1] >= 0.0)
            sum += "   ";
          else
            sum += "  ";
          sum += scientific_not(this->_array[i*_n + _n-1],3);
        }
        if (i < _m-1)
          sum += "\n  ";
      }
      sum += "    ...";
      if (_n < 10)
      {
        for (uint32_t j = 0; j < _n-1; j++)
        {
          sum += "         ...";
        }
        sum += "\n  ";
      }
      else
      {
        sum += "         ...         ...      ...   ...   ...       ...         ...         ...\n  ";
      }
      for (uint32_t i = _m-3; i < _m; i++)
      {
        if (_n < 10)
        {
          if (_array[i*_n] >= 0.0)
            sum += " ";
          for (uint32_t j = 0; j < _n; j++)
          {
            sum += scientific_not(this->_array[i*_n + j],3);
            if (j < _n-1)
            {
              if (_array[i*_n + j+1] >= 0.0)
                sum += "   ";
              else
                sum += "  ";
            }
          }
        }
        else
        {
          if (_array[i*_n] > 0.0)
            sum += " ";
          sum += scientific_not(this->_array[i*_n + 0],3);
          if (_array[i*_n + 1] >= 0.0)
            sum += "   ";
          else
            sum += "  ";
          sum += scientific_not(this->_array[i*_n + 1],3);
          if (_array[i*_n + 2] >= 0.0)
            sum += "   ";
          else
            sum += "  ";
          sum += scientific_not(this->_array[i*_n + 2],3);
          sum += "   ";
          sum += "...   ...   ...   ";
          if (_array[i*_n + _n-3] >= 0.0)
            sum += " ";
          sum += scientific_not(this->_array[i*_n + _n-3],3);
          if (_array[i*_n + _n-2] >= 0.0)
            sum += "   ";
          else
            sum += "  ";
          sum += scientific_not(this->_array[i*_n + _n-2],3);
          if (_array[i*_n + _n-1] >= 0.0)
            sum += "   ";
          else
            sum += "  ";
          sum += scientific_not(this->_array[i*_n + _n-1],3);
        }
        if (i < _m-1)
          sum += "\n  ";
      }
    }
    sum += "  ]";
    return sum;
  }
  template<typename T>
  Matrix<T> Matrix<T>::transpose() const
  {
    std::vector<T> new_array(_m*_n);
    for (uint32_t i = 0 ; i < _m; i++)
    {
      for (uint32_t j = 0; j < _n; j++)
      {
        new_array[j*_m + i] = _array[i*_n + j];
      }
    }
    std::string name = "(" + _name + ")^T";
    int m = _n;
    int n = _m;
    return Matrix<T>(name,m,n,new_array);
  }
  template<typename T>
  void Matrix<T>::transpose_inplace(bool inplace)
  {
    std::vector<T> new_array(_m*_n);
    for (uint32_t i = 0 ; i < _m; i++)
    {
      for (uint32_t j = 0; j < _n; j++)
      {
        new_array[j*_m + i] = _array[i*_n + j];
      }
    }
    std::string name = "(" + _name + ")^T";
    _name = name;
    int m = _n;
    int n = _m;
    _n = n;
    _m = m;
    _array = new_array;
  }
  template<typename T>
  T Matrix<T>::trace()
  {
    if (_m != _n)
    {
      std::cout << "Matrix is not square, trace is undefined!" << std::endl;
      return 0;
    }
    T result = 0;
    for (uint32_t i = 0; i < _m; i++)
    {
      result += _array[i*_n + i];
    }
    return result;
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Linear algebra tools
  //----------------------------------------------------------------------------
  template<typename T>
  bool Matrix<T>::isInvertible()
  {
  }
  template<typename T>
  void Matrix<T>::findSingularValues()
  {
  }
  template<typename T>
  Matrix<T> Matrix<T>::inverse()
  {
  }
  template<typename T>
  Matrix<T> Matrix<T>::pseudoInverse()
  {
  }
  template<typename T>
  std::tuple<Matrix<T>,Matrix<T>,Matrix<T>> Matrix<T>::LU()
  {
  }
  //  TODO: implement getL
  template<typename T>
  Matrix<T> Matrix<T>::getL(const Matrix<T>& perm)
  {
    return *this;
  }
  //  TODO: implement getU
  template<typename T>
  Matrix<T> Matrix<T>::getU(const Matrix<T>& perm)
  {
    return *this;
  }
  template<typename T>
  std::tuple<Matrix<T>,Matrix<T>> Matrix<T>::QR()
  {
  }
  template<typename T>
  std::tuple<Matrix<T>,Matrix<T>,Matrix<T>> Matrix<T>::SVD()
  {
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Various instantiators
  //----------------------------------------------------------------------------
  Matrix<double> identity_d(uint32_t m)
  {
    std::string name = "I_{" + std::to_string(m) + "x"
                     + std::to_string(m) + "}";
    Matrix<double> matrix(name, m, m, 0.0);

    for (uint32_t i = 0; i < m; i++)
    {
      matrix(i,i) = 1.0;
    }
    return matrix;
  }
  Matrix<double> zeros_d(uint32_t m)
  {
    Matrix<double> z(m, m, 0.0);
    return z;
  }
  Matrix<double> zeros_d(uint32_t m, uint32_t n)
  {
    Matrix<double> z(m, n, 0.0);
    return z;
  }
  Matrix<double> ones_d(uint32_t m)
  {
    Matrix<double> o(m, m, 1.0);
    return o;
  }
  Matrix<double> ones_d(uint32_t m, uint32_t n)
  {
    Matrix<double> o(m, n, 1.0);
    return o;
  }
  //----------------------------------------------------------------------------
  //  Permutation matrix
  //----------------------------------------------------------------------------
  Matrix<double> permutationMatrix_d(const uint32_t& m,
                              const std::vector<uint32_t> pivot)
  {
    //  generate a permutation matrix from a set of pivot indices
    std::vector<double> swaps(m*m, 0);
    std::vector<uint32_t> p(m,0);
    for (int i = 0; i < m; i++)
    {
      p[i] = i;
    }
    for (int i = 0; i < m; i++)
    {
      int temp;
      temp = p[pivot[i]-1];
      p[pivot[i]-1] = p[i];
      p[i] = temp;
    }
    for(uint32_t i = 0; i < m; i++)
    {
      swaps[p[i]*m + i] = 1;
    }
    std::string name = "(" + std::to_string(m) = "x"
                        + std::to_string(m) + ") perm";
    return Matrix<double>(name,m,m,swaps);
  }
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Overload of the ostream operator<<
  //  This has the same functionality as print, except it can be used in
  //  a std::cout statement.
  //----------------------------------------------------------------------------
  template<typename T>
  std::ostream& operator<<(std::ostream& os, const Matrix<T>& matrix)
  {
    os << "(" << matrix.getNumRows() << " x " << matrix.getNumCols() << ") Matrix";
    if (matrix.getName() != " ")
    {
      os << ": '" << matrix.getName() << "'";
    }
    os << "\n[ ";
    for (uint32_t i = 0; i < matrix.getNumRows(); i++)
    {
      os << "[ ";
      for (uint32_t j = 0; j < matrix.getNumCols(); j++)
      {
        os << matrix(i,j) << " ";
      }
      os << "]";
      if (i < matrix.getNumRows()-1)
      {
        os << "\n  ";
      }
    }
    os << " ]" << std::endl;
    return os;
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
  Vector<double> DGEMV(double& alpha, Matrix<double>& A,
                       Vector<double>& x)
  {
    //  container for the vector y
    std::vector<double> y(x.getDim());
    cblas_dgemv(CblasRowMajor, //  row major order
                CblasNoTrans,  //  don't transpose A
                A.getNumRows(),//  number of rows of A
                A.getNumCols(),//  number of columns of A
                alpha,         //  scaling factor for A*x
                A.data(),      //  pointer to A array
                A.getNumRows(),//  number of rows in A
                x.data(),      //  pointer to x vector
                1,             //  increment for elements of x
                1,             //  scaling factor for y
                y.data(),      //  pointer to y vector
                1);            //  increment for elements of y
    //  generate a new name for the product
    std::string name;
    if (alpha != 0)
    {
      name += std::to_string(alpha) + " * ";
    }
    name += A.getName() + " * " + x.getName();
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
  void DGEMV(double& alpha, Matrix<double>& A,
             Vector<double>& x, double& beta, Vector<double>& y)
  {
    cblas_dgemv(CblasRowMajor, //  row major order
                CblasNoTrans,  //  don't transpose A
                A.getNumRows(),//  number of rows of A
                A.getNumCols(),//  number of columns of A
                alpha,         //  scaling factor for A*x
                A.data(),      //  pointer to A array
                A.getNumRows(),//  number of rows in A
                x.data(),      //  pointer to x vector
                1,             //  increment for elements of x
                beta,          //  scaling factor for y
                y.data(),      //  pointer to y data
                1);            //  increment for the elements of y
    //  generate a new name for the product
    std::string name;
    if (alpha != 0)
    {
      name += std::to_string(alpha) + " * ";
    }
    name += A.getName() + " * " + x.getName();
    if (beta != 0)
    {
      name += " + " + std::to_string(beta) + " * " + y.getName();
    }
    y.setName(name);
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
  Matrix<double> DGER(double& alpha, Vector<double>& x, Vector<double>& y)
  {
    std::vector<double> A(x.getDim() * y.getDim());
    cblas_dger(CblasRowMajor,//  row major order
               x.getDim(),   //  number of rows of A
               y.getDim(),   //  number of columns of A
               alpha,        //  scaling factor for A*x
               x.data(),     //  pointer to x vector
               1,            //  increment for elements of x
               y.data(),     //  pointer to y data
               1,            //  increment for the elements of y
               A.data(),     //  pointer to the matrix A
               y.getDim());  //  leading dimension of A
    //  generate a new name for the product
    std::string name;
    if (alpha != 0 && x.getName() != " " && y.getName() != " ")
    {
      name += std::to_string(alpha) + " * " + x.getName()
              + " * " + y.getName() + "^T";
    }
    return Matrix<double>(name,x.getDim(),y.getDim(),A);
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
  void DGER(double& alpha, Vector<double>& x, Vector<double>& y,
                      Matrix<double>& m)
  {
    if (x.getDim() != m.getNumRows() || y.getDim() != m.getNumCols())
    {
      std::cout << "Vectors and Matrix are incompatible!" << std::endl;
      return;
    }
    cblas_dger(CblasRowMajor,//  row major order
               x.getDim(),   //  number of rows of A
               y.getDim(),   //  number of columns of A
               alpha,        //  scaling factor for A*x
               x.data(),     //  pointer to x vector
               1,            //  increment for elements of x
               y.data(),     //  pointer to y data
               1,            //  increment for the elements of y
               m.data(),     //  pointer to the matrix A
               y.getDim());  //  leading dimension of A
    //  generate a new name for the product
    std::string name;
    if (alpha != 0 && x.getName() != " " && y.getName() != " ")
    {
      name += std::to_string(alpha) + " * " + x.getName()
              + " * " + y.getName() + "^T";
    }
    if (m.getName() != " ")
    {
      name += " + " + m.getName();
    }
    m.setName(name);
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
  Matrix<double> DGEMM(const double& alpha, const Matrix<double>& A,
                       const Matrix<double>& B)
  {
    if (A.getNumCols() != B.getNumRows())
    {
      std::cout << "Matrices are incompatible!" << std::endl;
      return Matrix<double>("zeros",1,0.0);
    }
    std::vector<double> C(A.getNumRows() * B.getNumCols());
    cblas_dgemm(CblasRowMajor,  //  row major order
                CblasNoTrans,   //  don't tranpose A
                CblasNoTrans,   //  don't transpose B
                A.getNumRows(), //  number of rows of A
                B.getNumCols(), //  number of cols of B
                A.getNumCols(), //  number of cols of A
                alpha,          //  scaling factor for A*B
                A.data(),       //  pointer to elements of A
                A.getNumCols(), //  leading dimension of A
                B.data(),       //  pointer to elements of B
                B.getNumCols(), //  leading dimension of B
                1,              //  coefficient beta=1
                C.data(),       //  pointer to elements of C
                B.getNumCols());//  leading dimension of C
    //  generate a new name for the product
    std::string name;
    if (alpha != 0 && A.getName() != " " && B.getName() != " ")
    {
      name += std::to_string(alpha) + " * " + A.getName()
              + " * " + B.getName();
    }
    return Matrix<double>(name,A.getNumRows(),B.getNumCols(),C);
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
  void DGEMM(const double& alpha, const Matrix<double>& A,
             const Matrix<double>& B, const double& beta,
             Matrix<double>& C)
  {
    if (A.getNumCols() != B.getNumRows() || A.getNumRows() != C.getNumRows()
        || B.getNumCols() != C.getNumCols())
    {
      std::cout << "Matrices are incompatible!" << std::endl;
      return;
    }
    cblas_dgemm(CblasRowMajor,  //  row major order
                CblasNoTrans,   //  don't tranpose A
                CblasNoTrans,   //  don't transpose B
                A.getNumRows(), //  number of rows of A
                B.getNumCols(), //  number of cols of B
                A.getNumCols(), //  number of cols of A
                alpha,          //  scaling factor for A*B
                A.data(),       //  pointer to elements of A
                A.getNumCols(), //  leading dimension of A
                B.data(),       //  pointer to elements of B
                B.getNumCols(), //  leading dimension of B
                beta,           //  coefficient beta=1
                C.data(),       //  pointer to elements of C
                B.getNumCols());//  leading dimension of C
    //  generate a new name for the product
    std::string name;
    if (alpha != 0 && A.getName() != " " && B.getName() != " ")
    {
      name += std::to_string(alpha) + " * " + A.getName()
              + " * " + B.getName();
    }
    C.setName(name);
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DGEMM (no alpha) - generic matrix-matrix multiplication
  //  Arguments:  A     - (m x k)-matrix
  //              B     - (k x n)-matrix
  //
  //  Returns:    m * n  ( (m x n)-matrix )
  //----------------------------------------------------------------------------
  Matrix<double> DGEMM(const Matrix<double>& A,
                       const Matrix<double>& B)
  {
    if (A.getNumCols() != B.getNumRows())
    {
      std::cout << "Matrices are incompatible!" << std::endl;
      return Matrix<double>("zeros",1,0.0);
    }
    std::vector<double> C(A.getNumRows() * B.getNumCols());
    cblas_dgemm(CblasRowMajor,  //  row major order
                CblasNoTrans,   //  don't tranpose A
                CblasNoTrans,   //  don't transpose B
                A.getNumRows(), //  number of rows of A
                B.getNumCols(), //  number of cols of B
                A.getNumCols(), //  number of cols of A
                1.0,            //  scaling factor for A*B
                A.data(),       //  pointer to elements of A
                A.getNumCols(), //  leading dimension of A
                B.data(),       //  pointer to elements of B
                B.getNumCols(), //  leading dimension of B
                1.0,            //  coefficient beta=1
                C.data(),       //  pointer to elements of C
                B.getNumCols());//  leading dimension of C
    //  generate a new name for the product
    std::string name;
    if (A.getName() != " " && B.getName() != " ")
    {
      name += A.getName() + " * " + B.getName();
    }
    return Matrix<double>(name,A.getNumRows(),B.getNumCols(),C);
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  LAPACK methods
  //----------------------------------------------------------------------------
  //  DGELS - solve the linear system min_C|B - AC|
  //  Arguments:  A     - (m x k)-matrix
  //              B     - (m x n)-matrix
  //
  //  Returns:    C_min  ( (k x n)-matrix )
  //----------------------------------------------------------------------------
  Matrix<double> DGELS(const Matrix<double>& A, const Matrix<double>& B)
  {
    if (A.getNumRows() != B.getNumRows())
    {
      std::cout << "Matrices are incompatible!" << std::endl;
      return Matrix<double>("zeros",1,0.0);
    }
    Matrix<double> QR(A);
    std::vector<double> c = B.getArray();
    int info;
    info = LAPACKE_dgels(LAPACK_ROW_MAJOR,//  row major layout
                         'N',             //  don't transpose A
                         A.getNumRows(),  //  number of rows of A
                         A.getNumCols(),  //  number of columns of A
                         B.getNumCols(),  //  number of columns of B
                         QR.data(),       //  pointer to elements of A
                         A.getNumCols(),  //  leading dimension of A
                         c.data(),        //  pointer to elements of B
                         B.getNumCols()); //  leading dimension of B
    if( info > 0 )
    {
      std::cout << "The diagonal element " << info << " of the triangular "
                << "factor of A is zero, so that A does not have full rank;\n"
                << "the least squares solution could not be computed.\n";
    }
    //  Cut the result according to (A.getNumCols() x B.getNumCols())
    c.resize(A.getNumCols() * B.getNumCols());
    std::string name;
    if (A.getName() != " " && B.getName() != " ")
    {
      name += "min_C||" + A.getName() + "*C - " + B.getName() + "|";
    }
    else
    {
      name  = " ";
    }
    return Matrix<double>(name,A.getNumCols(),B.getNumCols(),c);
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DGELS - solve the linear system min_u|v - Au|
  //  Arguments:  A     - (m x k)-matrix
  //              v     - (m)-dim vector
  //
  //  Returns:    u_min  ( (k)-dim vector )
  //----------------------------------------------------------------------------
  Vector<double> DGELS(const Matrix<double>& A, const Vector<double>& v)
  {
    if (A.getNumRows() != v.getDim())
    {
      std::cout << "Matrices and vector are incompatible!" << std::endl;
      return Vector<double>("zeros",1,0.0);
    }
    Matrix<double> QR(A);
    //  ------------------------------------------------------------------------
    //  WARNING: while A can be a general m x n matrix, the vector u
    //  must have at least as many rows as the number of columns of A,
    //  otherwise we get a segfault from LAPACK.  To prevent this, u
    //  is initialized with _dim = A.getNumCols() with zeros.  Then, the
    //  first v.getDim() elements are filled with v's elements.
    //  (N. Carrara - 6/30/2020)
    //  Error fixed with commit - bc555791c2ebc9aadebe44e5b74fbc367b7e2123.
    //--------------------------------------------------------------------------
    std::vector<double> u(A.getNumCols(),0.0);
    for (uint32_t i = 0; i < v.getDim(); i++)
    {
      u[i] = v(i);
    }
    int info;
    info = LAPACKE_dgels(LAPACK_ROW_MAJOR,//  row major layout
                         'N',             //  don't transpose A
                         A.getNumRows(),  //  number of rows of A
                         A.getNumCols(),  //  number of rows of A
                         1,               //  dimension of v
                         QR.data(),       //  pointer to elements of A
                         A.getNumCols(),  //  leading dimension of A
                         u.data(),        //  pointer to elements of v
                         1);              //  leading dimension of v
    if( info > 0 )
    {
      std::cout << "The diagonal element " << info << " of the triangular "
                << "factor of A is zero, so that A does not have full rank;\n"
                << "the least squares solution could not be computed.\n";
    }
    //  Cut the result according to (A.getNumCols())
    u.resize(A.getNumCols());
    std::string name;
    if (A.getName() != " " && v.getName() != " ")
    {
      name += "min_u||" + A.getName() + "*u - " + v.getName() + "|";
    }
    else
    {
      name  = " ";
    }
    return Vector<double>(name,u);
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DGELSY - solve the linear system min_C|B - AC| using
  //           complete orthogonal factorization.
  //  Arguments:  A     - (m x k)-matrix
  //              B     - (m x n)-matrix
  //
  //  Returns:    C_min  ( (k x n)-matrix )
  //----------------------------------------------------------------------------
  Matrix<double> DGELSY(const Matrix<double>& A, const Matrix<double>& B)
  {
    if (A.getNumRows() != B.getNumRows())
    {
      std::cout << "Matrices are incompatible!" << std::endl;
      return Matrix<double>("zeros",1,0.0);
    }
    Matrix<double> QR(A);
    std::vector<double> c = B.getArray();
    int info;
    int jpvt[A.getNumCols()];
    double rcond;
    int rank;
    info = LAPACKE_dgelsy(LAPACK_ROW_MAJOR,//  row major layout
                          A.getNumRows(),  //  number of rows of A
                          A.getNumCols(),  //  number of columns of A
                          B.getNumCols(),  //  number of columns of B
                          QR.data(),       //  pointer to elements of A
                          A.getNumCols(),  //  leading dimension of A
                          c.data(),        //  pointer to elements of B
                          B.getNumCols(),  //  leading dimension of B
                          jpvt,            //  workspace array
                          rcond,           //  used for rank of A
                          &rank);          //  the effective rank of A
    if( info > 0 )
    {
      std::cout << "The diagonal element " << info << " of the triangular "
                << "factor of A is zero, so that A does not have full rank;\n"
                << "the least squares solution could not be computed.\n";
    }
    //  Cut the result according to (A.getNumCols() x B.getNumCols())
    c.resize(A.getNumCols() * B.getNumCols());
    std::string name;
    if (A.getName() != " " && B.getName() != " ")
    {
      name += "min_C||" + A.getName() + "*C - " + B.getName() + "|";
    }
    else
    {
      name  = " ";
    }
    //  set the effective rank of A
    A.setRank(rank);
    return Matrix<double>(name,A.getNumCols(),B.getNumCols(),c);
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DGELSY - solve the linear system min_u|v - Au|
  //           using complete orthogonal factorization
  //  Arguments:  A     - (m x k)-matrix
  //              v     - (m)-dim vector
  //
  //  Returns:    u_min  ( (k)-dim vector )
  //----------------------------------------------------------------------------
  Vector<double> DGELSY(const Matrix<double>& A, const Vector<double>& v)
  {
    if (A.getNumRows() != v.getDim())
    {
      std::cout << "Matrices and vector are incompatible!" << std::endl;
      return Vector<double>("zeros",1,0.0);
    }
    Matrix<double> QR(A);
    //  ------------------------------------------------------------------------
    //  WARNING: while A can be a general m x n matrix, the vector u
    //  must have at least as many rows as the number of columns of A,
    //  otherwise we get a segfault from LAPACK.  To prevent this, u
    //  is initialized with _dim = A.getNumCols() with zeros.  Then, the
    //  first v.getDim() elements are filled with v's elements.
    //  (N. Carrara - 6/30/2020)
    //  Error fixed with commit - bc555791c2ebc9aadebe44e5b74fbc367b7e2123.
    //--------------------------------------------------------------------------
    std::vector<double> u(A.getNumCols(),0.0);
    for (uint32_t i = 0; i < v.getDim(); i++)
    {
      u[i] = v(i);
    }    int info;
    int jpvt[A.getNumCols()];
    double rcond;
    int rank;
    info = LAPACKE_dgelsy(LAPACK_ROW_MAJOR,//  row major layout
                          A.getNumRows(),  //  number of rows of A
                          A.getNumCols(),  //  number of columns of A
                          1,               //  dimension of v
                          QR.data(),       //  pointer to elements of A
                          A.getNumCols(),  //  leading dimension of A
                          u.data(),        //  pointer to elements of v
                          1,               //  leading dimension of v
                          jpvt,            //  workspace array
                          rcond,           //  used for rank of A
                          &rank);          //  the effective rank of A
    if( info > 0 )
    {
      std::cout << "The diagonal element " << info << " of the triangular "
                << "factor of A is zero, so that A does not have full rank;\n"
                << "the least squares solution could not be computed.\n";
    }
    //  Cut the result according to (A.getNumCols())
    u.resize(A.getNumCols());
    std::string name;
    if (A.getName() != " " && v.getName() != " ")
    {
      name += "min_u||" + A.getName() + "*u - " + v.getName() + "|";
    }
    else
    {
      name  = " ";
    }
    //  set the effective rank of A
    A.setRank(rank);
    return Vector<double>(name,u);
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DGELSD - solve the linear system min_C|B - AC|
  //           using SVD and the divide and conquer method
  //  Arguments:  A     - (m x k)-matrix
  //              B     - (m x n)-matrix
  //
  //  Returns:    C_min  ( (k x n)-matrix )
  //----------------------------------------------------------------------------
  Matrix<double> DGELSD(const Matrix<double>& A, const Matrix<double>& B)
  {
    if (A.getNumRows() != B.getNumRows())
    {
      std::cout << "Matrices are incompatible!" << std::endl;
      return Matrix<double>("zeros",1,0.0);
    }
    Matrix<double> SVD(A);
    std::vector<double> c = B.getArray();
    int info;
    std::vector<double> singular(std::min(A.getNumRows(),A.getNumCols()));
    double rcond;
    int rank;
    info = LAPACKE_dgelsd(LAPACK_ROW_MAJOR,//  row major layout
                          A.getNumRows(),  //  number of rows of A
                          A.getNumCols(),  //  number of columns of A
                          B.getNumCols(),  //  number of columns of B
                          SVD.data(),      //  pointer to elements of A
                          A.getNumCols(),  //  leading dimension of A
                          c.data(),        //  pointer to elements of B
                          B.getNumCols(),  //  leading dimension of B
                          singular.data(), //  pointer to the singular values
                          rcond,           //  condition number
                          &rank);          //  effective rank of A.
    if( info > 0 )
    {
      std::cout << "The diagonal element " << info << " of the triangular "
                << "factor of A is zero, so that A does not have full rank;\n"
                << "the least squares solution could not be computed.\n";
    }
    //  Cut the result according to (A.getNumCols() x B.getNumCols())
    c.resize(A.getNumCols() * B.getNumCols());
    std::string name;
    if (A.getName() != " " && B.getName() != " ")
    {
      name += "min_C||" + A.getName() + "*C - " + B.getName() + "|";
    }
    else
    {
      name  = " ";
    }
    //  set the singular values of A
    A.setSingularValues(singular);
    //  set the effective rank of A
    A.setRank(rank);
    return Matrix<double>(name,A.getNumCols(),B.getNumCols(),c);
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DGELSD - solve the linear system min_u|v - Au|
  //           using SVD and the divide and conquer method
  //  Arguments:  A     - (m x k)-matrix
  //              v     - (m)-dim vector
  //
  //  Returns:    u_min  ( (k)-dim vector )
  //----------------------------------------------------------------------------
  Vector<double> DGELSD(const Matrix<double>& A, const Vector<double>& v)
  {
    if (A.getNumRows() != v.getDim())
    {
      std::cout << "Matrices and vector are incompatible!" << std::endl;
      return Vector<double>("zeros",1,0.0);
    }
    Matrix<double> SVD(A);
    //  ------------------------------------------------------------------------
    //  WARNING: while A can be a general m x n matrix, the vector u
    //  must have at least as many rows as the number of columns of A,
    //  otherwise we get a segfault from LAPACK.  To prevent this, u
    //  is initialized with _dim = A.getNumCols() with zeros.  Then, the
    //  first v.getDim() elements are filled with v's elements.
    //  (N. Carrara - 6/30/2020)
    //  Error fixed with commit - bc555791c2ebc9aadebe44e5b74fbc367b7e2123.
    //--------------------------------------------------------------------------
    std::vector<double> u(A.getNumCols(),0.0);
    for (uint32_t i = 0; i < v.getDim(); i++)
    {
      u[i] = v(i);
    }
    int info;
    std::vector<double> singular(std::min(A.getNumRows(),A.getNumCols()));
    double rcond;
    int rank;
    info = LAPACKE_dgelsd(LAPACK_ROW_MAJOR,//  row major layout
                          A.getNumRows(),  //  number of rows of A
                          A.getNumCols(),  //  number of columns of A
                          1,               //  dimension of v
                          SVD.data(),      //  pointer to elements of A
                          A.getNumCols(),  //  leading dimension of A
                          u.data(),        //  pointer to elements of v
                          1,               //  leading dimension of v
                          singular.data(), //  pointer to the singular values
                          rcond,           //  condition number
                          &rank);          //  effective rank of A.
    if( info > 0 )
    {
      std::cout << "The diagonal element " << info << " of the triangular "
                << "factor of A is zero, so that A does not have full rank;\n"
                << "the least squares solution could not be computed.\n";
    }
    //  Cut the result according to (A.getNumCols())
    u.resize(A.getNumCols());
    std::string name;
    if (A.getName() != " " && v.getName() != " ")
    {
      name += "min_u||" + A.getName() + "*u - " + v.getName() + "|";
    }
    else
    {
      name  = " ";
    }
    //  set the singular values of A
    A.setSingularValues(singular);
    //  set the effective rank of A
    A.setRank(rank);
    return Vector<double>(name,u);
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DGELSS - solve the linear system min_C|B - AC|
  //           using SVD.
  //  Arguments:  A     - (m x k)-matrix
  //              B     - (m x n)-matrix
  //
  //  Returns:    C_min  ( (k x n)-matrix )
  //----------------------------------------------------------------------------
  Matrix<double> DGELSS(const Matrix<double>& A, const Matrix<double>& B)
  {
    if (A.getNumRows() != B.getNumRows())
    {
      std::cout << "Matrices are incompatible!" << std::endl;
      return Matrix<double>("zeros",1,0.0);
    }
    Matrix<double> SVD(A);
    std::vector<double> c = B.getArray();
    int info;
    std::vector<double> singular(std::min(A.getNumRows(),A.getNumCols()));
    double rcond;
    int rank;
    info = LAPACKE_dgelss(LAPACK_ROW_MAJOR,//  row major layout
                          A.getNumRows(),  //  number of rows of A
                          A.getNumCols(),  //  number of columns of A
                          B.getNumCols(),  //  number of columns of B
                          SVD.data(),      //  pointer to elements of A
                          A.getNumCols(),  //  leading dimension of A
                          c.data(),        //  pointer to elements of B
                          B.getNumCols(),  //  leading dimension of B
                          singular.data(), //  pointer to the singular values
                          rcond,           //  condition number
                          &rank);          //  effective rank of A.
    if( info > 0 )
    {
      std::cout << "The diagonal element " << info << " of the triangular "
                << "factor of A is zero, so that A does not have full rank;\n"
                << "the least squares solution could not be computed.\n";
    }
    //  Cut the result according to (A.getNumCols() x B.getNumCols())
    c.resize(A.getNumCols() * B.getNumCols());
    std::string name;
    if (A.getName() != " " && B.getName() != " ")
    {
      name += "min_C||" + A.getName() + "*C - " + B.getName() + "|";
    }
    else
    {
      name  = " ";
    }
    //  set the singular values of A
    A.setSingularValues(singular);
    //  set the effective rank of A
    A.setRank(rank);
    return Matrix<double>(name,A.getNumCols(),B.getNumCols(),c);
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DGELSS - solve the linear system min_u|v - Au|
  //           using SVD.
  //  Arguments:  A     - (m x k)-matrix
  //              v     - (m)-dim vector
  //
  //  Returns:    u_min  ( (k)-dim vector )
  //----------------------------------------------------------------------------
  Vector<double> DGELSS(const Matrix<double>& A, const Vector<double>& v)
  {
    if (A.getNumRows() != v.getDim())
    {
      std::cout << "Matrices and vector are incompatible!" << std::endl;
      return Vector<double>("zeros",1,0.0);
    }
    Matrix<double> SVD(A);
    //  ------------------------------------------------------------------------
    //  WARNING: while A can be a general m x n matrix, the vector u
    //  must have at least as many rows as the number of columns of A,
    //  otherwise we get a segfault from LAPACK.  To prevent this, u
    //  is initialized with _dim = A.getNumCols() with zeros.  Then, the
    //  first v.getDim() elements are filled with v's elements.
    //  (N. Carrara - 6/30/2020)
    //  Error fixed with commit - bc555791c2ebc9aadebe44e5b74fbc367b7e2123.
    //--------------------------------------------------------------------------
    std::vector<double> u(A.getNumCols(),0.0);
    for (uint32_t i = 0; i < v.getDim(); i++)
    {
      u[i] = v(i);
    }
    int info;
    std::vector<double> singular(std::min(A.getNumRows(),A.getNumCols()));
    double rcond;
    int rank;
    info = LAPACKE_dgelss(LAPACK_ROW_MAJOR,//  row major layout
                          A.getNumRows(),  //  number of rows of A
                          A.getNumCols(),  //  number of columns of A
                          1,               //  dimension of v
                          SVD.data(),      //  pointer to elements of A
                          A.getNumCols(),  //  leading dimension of A
                          u.data(),        //  pointer to elements of v
                          1,               //  leading dimension of v
                          singular.data(), //  pointer to the singular values
                          rcond,           //  condition number
                          &rank);          //  effective rank of A.
    if( info > 0 )
    {
      std::cout << "The diagonal element " << info << " of the triangular "
                << "factor of A is zero, so that A does not have full rank;\n"
                << "the least squares solution could not be computed.\n";
    }
    //  Cut the result according to (A.getNumCols())
    u.resize(A.getNumCols());
    std::string name;
    if (A.getName() != " " && v.getName() != " ")
    {
      name += "min_u||" + A.getName() + "*u - " + v.getName() + "|";
    }
    else
    {
      name  = " ";
    }
    //  set the singular values of A
    A.setSingularValues(singular);
    //  set the effective rank of A
    A.setRank(rank);
    return Vector<double>(name,u);
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DGETRF - construct an P^-1LU factorization
  //  Arguments:  A     - (m x m)-matrix
  //
  //  Returns:    std::vector<uint32_t> (list of pivots for P^-1)
  //----------------------------------------------------------------------------
  std::vector<uint32_t> DGETRF(const Matrix<double>& A)
  {
    if (A.getNumRows() != A.getNumCols())
    {
      std::cout << "Matrix is not square!" << std::endl;
      return std::vector<uint32_t>(0);
    }
    Matrix<double> LU(A);
    std::vector<uint32_t> ipiv(std::min(A.getNumRows(),A.getNumCols()));
    int info;
    info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR,// row major layout
                          A.getNumRows(),  // number of rows of A
                          A.getNumCols(),  // number of columns of A
                          LU.data(),       // pointer to the elements of A
                          A.getNumCols(),  // leading dimension of A
                          ipiv.data());    // pointer to pivot list
    return ipiv;
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
  //                    [l_m1  l_m2  ---  u_mm]]
  //----------------------------------------------------------------------------
  Matrix<double> DGETRF_L_U(const Matrix<double>& A)
  {
    if (A.getNumRows() != A.getNumCols())
    {
      std::cout << "Matrix is not square!" << std::endl;
      return Matrix<double>("zeros",1,0.0);
    }
    Matrix<double> LU(A);
    std::vector<uint32_t> ipiv(std::min(A.getNumRows(),A.getNumCols()));
    int info;
    info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR,// row major layout
                          A.getNumRows(),  // number of rows of A
                          A.getNumCols(),  // number of columns of A
                          LU.data(),       // pointer to the elements of A
                          A.getNumCols(),  // leading dimension of A
                          ipiv.data());    // pointer to pivot list
    std::string name = "LU-Matrix";
    if (A.getName() != " ")
    {
      name += " of " + A.getName();
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
  std::tuple<Matrix<double>,Matrix<double>> DGETRF_LU(const Matrix<double>& A)
  {
    if (A.getNumRows() != A.getNumCols())
    {
      std::cout << "Matrix is not square!" << std::endl;
      return {Matrix<double>("zeros",1,0.0),Matrix<double>("zeros",1,0.0)};
    }
    Matrix<double> LU(A);
    std::vector<uint32_t> ipiv(std::min(A.getNumRows(),A.getNumCols()));
    int info;
    info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR,// row major layout
                          A.getNumRows(),  // number of rows of A
                          A.getNumCols(),  // number of columns of A
                          LU.data(),       // pointer to the elements of A
                          A.getNumCols(),  // leading dimension of A
                          ipiv.data());    // pointer to pivot list
    std::vector<double> l(A.getNumRows()*A.getNumRows(),0.0);
    std::vector<double> u(A.getNumRows()*A.getNumRows(),0.0);
    for (uint32_t i = 0; i < A.getNumRows(); i++)
    {
      for (uint32_t j = 0; j < A.getNumRows(); j++)
      {
        if (i == j)
        {
          l[i*A.getNumRows() + j] = 1.0;
          u[i*A.getNumRows() + j] = LU(i,j);
        }
        else if (j < i)
        {
          l[i*A.getNumRows() + j] = LU(i,j);
        }
        else
        {
          u[i*A.getNumRows() + j] = LU(i,j);
        }
      }
    }
    std::string name_l = "L-Matrix";
    std::string name_u = "U-Matrix";
    if (A.getName() != " ")
    {
      name_l += " of " + A.getName();
      name_u += " of " + A.getName();
    }
    Matrix<double> L(name_l,A.getNumRows(),l);
    Matrix<double> U(name_u,A.getNumRows(),u);
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
  DGETRF_PLU(const Matrix<double>& A)
  {
    if (A.getNumRows() != A.getNumCols())
    {
      std::cout << "Matrix is not square!" << std::endl;
      return {Matrix<double>("zeros",1,0.0),
              Matrix<double>("zeros",1,0.0),
              Matrix<double>("zeros",1,0.0)};
    }
    Matrix<double> LU(A);
    std::vector<uint32_t> ipiv(std::min(A.getNumRows(),A.getNumCols()));
    int info;
    info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR,// row major layout
                          A.getNumRows(),  // number of rows of A
                          A.getNumCols(),  // number of columns of A
                          LU.data(),       // pointer to the elements of A
                          A.getNumCols(),  // leading dimension of A
                          ipiv.data());    // pointer to pivot list
    std::vector<double> l(A.getNumRows()*A.getNumRows(),0.0);
    std::vector<double> u(A.getNumRows()*A.getNumRows(),0.0);
    for (uint32_t i = 0; i < A.getNumRows(); i++)
    {
      for (uint32_t j = 0; j < A.getNumRows(); j++)
      {
        if (i == j)
        {
          l[i*A.getNumRows() + j] = 1.0;
          u[i*A.getNumRows() + j] = LU(i,j);
        }
        else if (j < i)
        {
          l[i*A.getNumRows() + j] = LU(i,j);
        }
        else
        {
          u[i*A.getNumRows() + j] = LU(i,j);
        }
      }
    }
    std::string name_p = "P-Matrix";
    std::string name_l = "L-Matrix";
    std::string name_u = "U-Matrix";
    if (A.getName() != " ")
    {
      name_p += " of " + A.getName();
      name_l += " of " + A.getName();
      name_u += " of " + A.getName();
    }
    Matrix<double> P = permutationMatrix_d(A.getNumRows(),ipiv);
    P.setName(name_p);
    Matrix<double> L(name_l,A.getNumRows(),l);
    Matrix<double> U(name_u,A.getNumRows(),u);
    return {P,L,U};
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  DGEQRF - compute a QR decomposition
  //  Arguments:  A     - (m x n)-matrix
  //
  //  Returns:    Vector<double> (elementary reflectors)
  //----------------------------------------------------------------------------
  Vector<double> DGEQRF(const Matrix<double>& A)
  {
    Matrix<double> QR(A);
    std::vector<double> reflectors(std::min(A.getNumRows(),A.getNumCols()));
    int info;
    info = LAPACKE_dgeqrf(LAPACK_ROW_MAJOR,  // row major order
                          A.getNumRows(),    // number of rows of A
                          A.getNumCols(),    // number of columns of A
                          QR.data(),         // pointer to the elements of A
                          A.getNumCols(),    // leading dimension of A
                          reflectors.data());// pointer to reflectors
    std::string name = "QR-reflectors";
    if (A.getName() != " ")
    {
      name += " of " + A.getName();
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
  Matrix<double> DORGQR(const Matrix<double>& A,
                             const Vector<double>& ref)
  {
    Matrix<double> QR(A);
    int info;
    info = LAPACKE_dorgqr(LAPACK_ROW_MAJOR,// row major order
                          A.getNumRows(),  // number of rows of A
                          A.getNumCols(),  // number of columns of A
                          ref.getDim(),    // number of elementary reflectors
                          QR.data(),       // pointer to the elements of A
                          A.getNumCols(),  // leading dimension of A
                          ref.data());     // pointer to reflectors
    std::string name = "Q-Matrix";
    if (A.getName() != " ")
    {
      name += " of " + A.getName();
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
  DGEQRF_QR(const Matrix<double>& A)
  {
    Matrix<double> Q(A);
    std::vector<double> reflectors(std::min(A.getNumRows(),A.getNumCols()));
    int info;
    info = LAPACKE_dgeqrf(LAPACK_ROW_MAJOR,  // row major order
                          A.getNumRows(),    // number of rows of A
                          A.getNumCols(),    // number of columns of A
                          Q.data(),          // pointer to the elements of A
                          A.getNumCols(),    // leading dimension of A
                          reflectors.data());// pointer to reflectors
    //  find Q using dorgqr
    info = LAPACKE_dorgqr(LAPACK_ROW_MAJOR,  // row major order
                          A.getNumRows(),    // number of rows of A
                          A.getNumCols(),    // number of columns of A
                          reflectors.size(), // number of elementary reflectors
                          Q.data(),          // pointer to the elements of A
                          A.getNumCols(),    // leading dimension of A
                          reflectors.data());// pointer to reflectors
    std::string name_q = "Q-Matrix";
    std::string name_r = "R-Matrix";
    if (A.getName() != " ")
    {
      name_q += " of " + A.getName();
      name_r += " of " + A.getName();
    }
    Q.setName(name_q);
    Matrix<double> Q_T(Q);
    Q_T.transpose_inplace();
    //  find R = Q^T*A
    Matrix<double> R = Q_T * A;
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
  Vector<double> DGESVD(const Matrix<double>& A)
  {
    Matrix<double> A_copy(A);
    Matrix<double> U(A.getNumRows());
    Matrix<double> VT(A.getNumCols());
    double superb[std::min(A.getNumRows(),A.getNumCols())-1];
    std::vector<double> singular(std::min(A.getNumRows(),A.getNumCols()));
    int info;
    info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR,// row major format
                          'N',             // no columns returned in U
                          'N',             // no rows returned in VT
                          A.getNumRows(),  // number of rows of A
                          A.getNumCols(),  // number of columns of A
                          A_copy.data(),   // pointer to the elements of A
                          A.getNumCols(),  // leading dimension of A
                          singular.data(), // pointer to singular values
                          U.data(),        // pointer to elements of U
                          U.getNumCols(),  // leading dimension of U
                          VT.data(),       // pointer to elements of VT
                          VT.getNumCols(), // leading dimension of VT
                          superb);         // workspace for subroutines
    std::string name = "Singular values";
    if (A.getName() != " ")
    {
      name += " of " + A.getName();
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
  DGESVD_SVD(const Matrix<double>& A)
  {
    Matrix<double> A_copy(A);
    Matrix<double> U(A.getNumRows());
    Matrix<double> VT(A.getNumCols());
    double superb[std::min(A.getNumRows(),A.getNumCols())-1];
    std::vector<double> singular(std::min(A.getNumRows(),A.getNumCols()));
    int info;
    info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR,// row major format
                          'A',             // all m rows returned in U
                          'A',             // all n rows returned in VT
                          A.getNumRows(),  // number of rows of A
                          A.getNumCols(),  // number of columns of A
                          A_copy.data(),   // pointer to the elements of A
                          A.getNumCols(),  // leading dimension of A
                          singular.data(), // pointer to singular values
                          U.data(),        // pointer to elements of U
                          U.getNumCols(),  // leading dimension of U
                          VT.data(),       // pointer to elements of VT
                          VT.getNumCols(), // leading dimension of VT
                          superb);         // workspace for subroutines
    std::string name_u = "U-Matrix";
    std::string name_s = "Sigma-Matrix";
    std::string name_vt = "(V)^T-Matrix";
    if (A.getName() != " ")
    {
      name_u += " of " + A.getName();
      name_s += " of " + A.getName();
      name_vt += " of " + A.getName();
    }
    U.setName(name_u);
    VT.setName(name_vt);
    Matrix<double> Sigma(name_s,A.getNumRows(),A.getNumCols(),0.0);
    for (uint32_t i = 0; i < singular.size(); i++)
    {
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
  Vector<double> DGESDD(const Matrix<double>& A)
  {
    Matrix<double> A_copy(A);
    Matrix<double> U(A.getNumRows());
    Matrix<double> VT(A.getNumCols());
    double superb[std::min(A.getNumRows(),A.getNumCols())-1];
    std::vector<double> singular(std::min(A.getNumRows(),A.getNumCols()));
    int info;
    info = LAPACKE_dgesdd(LAPACK_ROW_MAJOR,// row major format
                          'N',             // no columns returned in U or V
                          A.getNumRows(),  // number of rows of A
                          A.getNumCols(),  // number of columns of A
                          A_copy.data(),   // pointer to the elements of A
                          A.getNumCols(),  // leading dimension of A
                          singular.data(), // pointer to singular values
                          U.data(),        // pointer to elements of U
                          U.getNumCols(),  // leading dimension of U
                          VT.data(),       // pointer to elements of VT
                          VT.getNumCols());// leading dimension of VT
    std::string name = "Singular values";
    if (A.getName() != " ")
    {
      name += " of " + A.getName();
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
  DGESDD_SVD(const Matrix<double>& A)
  {
    Matrix<double> A_copy(A);
    Matrix<double> U(A.getNumRows());
    Matrix<double> VT(A.getNumCols());
    double superb[std::min(A.getNumRows(),A.getNumCols())-1];
    std::vector<double> singular(std::min(A.getNumRows(),A.getNumCols()));
    int info;
    info = LAPACKE_dgesdd(LAPACK_ROW_MAJOR,// row major format
                          'A',             // all m rows returned in U and VT
                          A.getNumRows(),  // number of rows of A
                          A.getNumCols(),  // number of columns of A
                          A_copy.data(),   // pointer to the elements of A
                          A.getNumCols(),  // leading dimension of A
                          singular.data(), // pointer to singular values
                          U.data(),        // pointer to elements of U
                          U.getNumCols(),  // leading dimension of U
                          VT.data(),       // pointer to elements of VT
                          VT.getNumCols());// leading dimension of VT
    std::string name_u = "U-Matrix";
    std::string name_s = "Sigma-Matrix";
    std::string name_vt = "(V)^T-Matrix";
    if (A.getName() != " ")
    {
      name_u += " of " + A.getName();
      name_s += " of " + A.getName();
      name_vt += " of " + A.getName();
    }
    U.setName(name_u);
    VT.setName(name_vt);
    Matrix<double> Sigma(name_s,A.getNumRows(),A.getNumCols(),0.0);
    for (uint32_t i = 0; i < singular.size(); i++)
    {
      Sigma(i,i) = singular[i];
    }
    return {U,Sigma,VT};
  }
  //----------------------------------------------------------------------------
}
