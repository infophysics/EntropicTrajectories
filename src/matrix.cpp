//  Matrix wrapper cpp file
#include "vector.h"
#include "matrix.h"

namespace ET
{

  template<typename T>
  Matrix<T>::Matrix() : _name(" ") {}

  template<typename T>
  Matrix<T>::~Matrix() {}

  //  Copy constructor
  template<typename T>
  Matrix<T>::Matrix(const Matrix<T>& m)
  {
    _mat = m._mat;
    _n = m.get_rows();
    _m = m.get_cols();
    _name = m.get_name();
  }

  template<typename T>
  Matrix<T>::Matrix(unsigned int n) : _n(n), _m(n), _name(" ")
  {
    _mat.resize(_n*_n,0.0);
  }

  template<typename T>
  Matrix<T>::Matrix(std::string name, unsigned int n)
  : _n(n), _m(n), _name(name)
  {
    _mat.resize(_n*_n,0.0);
  }

  template<typename T>
  Matrix<T>::Matrix(unsigned int n, unsigned int m) : _n(n), _m(m), _name(" ")
  {
    _mat.resize(_n*_m,0.0);
  }

  template<typename T>
  Matrix<T>::Matrix(std::string name, unsigned int n, unsigned int m)
  : _n(n), _m(m), _name(name)
  {
    _mat.resize(_n*_m,0.0);
  }

  template<typename T>
  Matrix<T>::Matrix(unsigned int n, unsigned int m, const T& init)
  : _n(n), _m(m), _name(" ")
  {
    _mat.resize(_n*_m, init);
  }

  template<typename T>
  Matrix<T>::Matrix(std::string name, unsigned int n,
    unsigned int m, const T& init)
  : _n(n), _m(m), _name(name)
  {
    _mat.resize(_n*_m, init);
  }

  template<typename T>
  Matrix<T>::Matrix(unsigned int n, std::vector<T> flat)
  : _n(n), _m(n), _name(" "), _mat(flat)
  {

  }
  template<typename T>
  Matrix<T>::Matrix(std::string name, unsigned int n, std::vector<T> flat)
  : _n(n), _m(n), _name(name), _mat(flat)
  {

  }
  template<typename T>
  Matrix<T>::Matrix(unsigned int n, unsigned int m, std::vector<T> flat)
  : _n(n), _m(m), _name(" "), _mat(flat)
  {

  }

  template<typename T>
  Matrix<T>::Matrix(std::string name, unsigned int n, unsigned int m, std::vector<T> flat)
  : _n(n), _m(m), _name(name), _mat(flat)
  {
  }
  template<typename T>
  unsigned int Matrix<T>::get_rows() const
  {
    return _n;
  }

  template<typename T>
  unsigned int Matrix<T>::get_cols() const
  {
    return _m;
  }

  template<typename T>
  std::string Matrix<T>::get_name() const
  {
    return _name;
  }

  template<typename T>
  std::vector<T> Matrix<T>::get_mat() const
  {
    return _mat;
  }

  template<typename T>
  std::vector<T> Matrix<T>::get_row(unsigned int i)
  {
    std::vector<T> row(_m,0.0);
    for (unsigned int j = 0; j < _m; j++) {
      row[j] = _mat[i*_m + j];
    }
    return row;
  }

  template<typename T>
  std::vector<T> Matrix<T>::get_col(unsigned int i)
  {
    std::vector<T> col(_n,0.0);
    for (unsigned int j = 0; j < _n; j++) {
      col[j] = _mat[j*_m + i];
    }
    return col;
  }

  //  Setters
  template<typename T>
  void Matrix<T>::set_name(std::string name)
  {
    _name = name;
  }
  template<typename T>
  void Matrix<T>::set_row(unsigned int i, std::vector<T> row)
  {
    for (unsigned int j = 0; j < _m; j++) {
      _mat[i*_m + j] = row[j];
    }
  }
  template<typename T>
  void Matrix<T>::set_col(unsigned int i, std::vector<T> col)
  {
    for (unsigned int j = 0; j < _n; j++) {
      _mat[j*_m + i] = col[j];
    }
  }

  //  Operator overloads
  template<typename T>
  Matrix<T>& Matrix<T>::operator=(const Matrix<T>& m)
  {
    if (&m == this)
      return *this;

    _n = m.get_rows();
    _m = m.get_cols();
    _mat.resize(_n*_m);
    for (unsigned int i = 0; i < _n*_m; i++) {
        _mat[i] = m(i);
    }
    return *this;
  }
  template<typename T>
  bool Matrix<T>::operator==(const Matrix<T>& m)
  {
    if (_n != m.get_rows() || _m != m.get_cols())
      return false;
    for (unsigned int i = 0; i < _n*_m; i++) {
        if (m(i) != _mat[i])
          return false;
    }
    return true;
  }
  //  Matrix algebra
  template<typename T>
  Matrix<T> Matrix<T>::operator+(const Matrix<T>& m)
  {
    if(_n != m.get_rows() || _m != m.get_cols())
      return *this;
    Matrix<T> l(_n, _m, 0.0);
    for (unsigned int i = 0; i < _n*_m; i++) {
        l(i) = _mat[i] + m(i);
    }
    return l;
  }
  template<typename T>
  Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& m)
  {
    if(_n != m.get_rows() || _m != m.get_cols())
      return *this;
    for (unsigned int i = 0; i < _n*_m; i++) {
        _mat[i] += m(i);
    }
    return *this;
  }
  template<typename T>
  Matrix<T> Matrix<T>::operator-(const Matrix<T>& m)
  {
    if(_n != m.get_rows() || _m != m.get_cols())
      return *this;
    Matrix<T> l(_n, _m, 0.0);
    for (unsigned int i = 0; i < _n*_m; i++) {
        l(i) = _mat[i] - m(i);
    }
    return l;
  }
  template<typename T>
  Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& m)
  {
    if(_n != m.get_rows() || _m != m.get_cols())
      return *this;
    for (unsigned int i = 0; i < _n*_m; i++) {
        _mat[i] -= m(i);
    }
    return *this;
  }
  template<typename T>
  Matrix<T> Matrix<T>::operator*(const Matrix<T>& m)
  {
    if(_m != m.get_rows())
      return *this;
    std::vector<T> l(_n*m.get_cols(),0.0);
    //  CBLAS function for matrix multiplication, A*B = C.
    /*  clbas_dgemm(Order  - either CblasRowMajor or CblasColumnMajor
                    TransA - either transpose or not for Matrix A
                    TransB - either transpose or not for Matrix B
                    M      - number of rows in A and C
                    N      - number of columns in B and C
                    K      - number of columns in A, number of rows in B
                    alpha  - scaling factor for A and B
                    A      - Matrix A, i.e. ref. to the beginning of the vector
                    Ida    - Size of first dimension of A; _n
                    B      - Matrix B, ref. to the beginning of B
                    Idb    - Size of the first dimension of B; m.get_cols()
                    beta   - scaling factor for C
                    C      - Matrix C, ref. to beginning of C
                    Idc    - Size of first dimension of C; _n
        Since we are using std::vector, we must pass in a reference to
        the pointer given by .begin().
    */
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                _n, m.get_cols(), _m, 1.0,
                & *_mat.begin(), _n,
                & *m.get_mat().begin(), _m,
                0.0, & *l.begin(), _n);
    return Matrix<T>(_n,m.get_cols(),l);
  }
  template<typename T>
  Matrix<T> Matrix<T>::brute_mul(const Matrix<T>& m)
  {
    if(_m != m.get_rows())
      return *this;
    Matrix<T> l(_n,m.get_cols(),0.0);
    for (unsigned int i = 0; i < _n; i++) {
      for (unsigned int j = 0; j < _m; j++) {
        for (unsigned int k = 0; k < _n; k++) {
          l(i,j) += this->_mat[i*_m + k] * m(k,j);
        }
      }
    }
    return l;
  }
  template<typename T>
  Matrix<T>& Matrix<T>::operator*=(const Matrix<T>& m)
  {
    Matrix<T> l = (*this) * m;
    (*this) = l;
    return *this;
  }
  //  Scalar algebra
  template<typename T>
  Matrix<T> Matrix<T>::operator+(const T& s)
  {
    Matrix<T> l(_n,_m,0.0);
    for (unsigned int i = 0; i < _n*_m; i++) {
        l(i) = _mat[i] + s;
    }
    return l;
  }
  template<typename T>
  Matrix<T> Matrix<T>::operator-(const T& s)
  {
    Matrix<T> l(_n,_m,0.0);
    for (unsigned int i = 0; i < _n*_m; i++) {
        l(i) = _mat[i] - s;
    }
    return l;
  }
  template<typename T>
  Matrix<T> Matrix<T>::operator*(const T& s)
  {
    Matrix<T> l(_n,_m,0.0);
    for (unsigned int i = 0; i < _n*_m; i++) {
        l(i) = _mat[i] * s;
    }
    return l;
  }
  template<typename T>
  Matrix<T> Matrix<T>::operator/(const T& s)
  {
    Matrix<T> l(_n,_m,0.0);
    if (s == 0)
      return l;
    for (unsigned int i = 0; i < _n*_m; i++) {
        l(i) = _mat[i] + s;
    }
    return l;
  }
  template<typename T>
  Matrix<T>& Matrix<T>::operator+=(const T& s)
  {
    for (unsigned int i = 0; i < _n*_m; i++) {
        _mat[i] += s;
    }
    return *this;
  }
  template<typename T>
  Matrix<T>& Matrix<T>::operator-=(const T& s)
  {
    for (unsigned int i = 0; i < _n*_m; i++) {
        _mat[i] -= s;
    }
    return *this;
  }
  template<typename T>
  Matrix<T>& Matrix<T>::operator*=(const T& s)
  {
    for (unsigned int i = 0; i < _n*_m; i++) {
        _mat[i] *= s;
    }
    return *this;
  }
  template<typename T>
  Matrix<T>& Matrix<T>::operator/=(const T& s)
  {
    if (s == 0)
      return *this;
    for (unsigned int i = 0; i < _n*_m; i++) {
        _mat[i] += s;
    }
    return *this;
  }
  //  Multiplying a vector
  template<typename T>
  Vector<T> Matrix<T>::operator*(const Vector<T>& v)
  {
    std::vector<T> vec(_n,0.0);
    Vector<T> v2(vec);
    if (_m != v.get_dim())
      return v2;
    for (unsigned int i = 0; i < _n; i++) {
      T temp = 0.0;
      for (unsigned int j = 0; j < _m; j++) {
        temp += this->_mat[i*_m + j] * v(j);
      }
      v2(i) = temp;
    }
    return v2;
  }
  //  Access Matrix<T>::operators
  template<typename T>
  T& Matrix<T>::operator()(const unsigned int& i, const unsigned int& j)
  {
    return this->_mat[i*_m + j];
  }
  template<typename T>
  const T& Matrix<T>::operator()(const unsigned int& i, const unsigned int& j) const
  {
    return this->_mat[i*_m + j];
  }
  template<typename T>
  T& Matrix<T>::operator()(const unsigned int& i)
  {
    return this->_mat[i];
  }
  template<typename T>
  const T& Matrix<T>::operator()(const unsigned int& i) const
  {
    return this->_mat[i];
  }
  template<typename T>
  std::ostream& operator<<(std::ostream& os, const Matrix<T>& m)
  {
    os << "(" << m.get_rows() << " x " << m.get_cols() << ") Matrix";
    if (m.get_name() != " ")
      os << ": '" << m.get_name() << "'";
    os << "\n[ ";
    for (unsigned int i = 0; i < m.get_rows(); i++) {
      os << "[ ";
      for (unsigned int j = 0; j < m.get_cols(); j++) {
        os << m(i,j) << " ";
      }
      os << "]";
      if (i < m.get_rows()-1)
        os << "\n  ";
    }
    os << " ]" << std::endl;
    return os;
  }

  //  Various methods
  template<typename T>
  void Matrix<T>::print()
  {
    std::cout << "(" << _n << " x " << _m << ") Matrix";
    if (_name != " ")
      std::cout << ": '" << _name << "'";
    std:: cout << "\n[ ";
    for (unsigned int i = 0; i < _n; i++) {
      std::cout << "[ ";
      for (unsigned int j = 0; j < _m; j++) {
        std::cout << this->_mat[i*_m + j] << " ";
      }
      std::cout << "]";
      if (i < _n-1)
        std::cout << "\n  ";
    }
    std::cout << " ]" << std::endl;
  }
  template<typename T>
  Matrix<T> Matrix<T>::transpose()
  {

  }
}
