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
    _mat.resize(_n);
    for (unsigned int i = 0; i < _mat.size(); i++) {
      _mat[i].resize(_m, 0.0);
    }
  }

  template<typename T>
  Matrix<T>::Matrix(std::string name, unsigned int n)
  : _n(n), _m(n), _name(name)
  {
    _mat.resize(_n);
    for (unsigned int i = 0; i < _mat.size(); i++) {
      _mat[i].resize(_m, 0.0);
    }
  }

  template<typename T>
  Matrix<T>::Matrix(unsigned int n, unsigned int m) : _n(n), _m(m), _name(" ")
  {
    _mat.resize(_n);
    for (unsigned int i = 0; i < _mat.size(); i++) {
      _mat[i].resize(_m, 0.0);
    }
  }

  template<typename T>
  Matrix<T>::Matrix(std::string name, unsigned int n, unsigned int m)
  : _n(n), _m(m), _name(name)
  {
    _mat.resize(_n);
    for (unsigned int i = 0; i < _mat.size(); i++) {
      _mat[i].resize(_m, 0.0);
    }
  }

  template<typename T>
  Matrix<T>::Matrix(unsigned int n, unsigned int m, const T& init)
  : _n(n), _m(m), _name(" ")
  {
    _mat.resize(_n);
    for (unsigned int i = 0; i < _mat.size(); i++) {
      _mat[i].resize(_m, init);
    }
  }

  template<typename T>
  Matrix<T>::Matrix(std::string name, unsigned int n,
    unsigned int m, const T& init)
  : _n(n), _m(m), _name(name)
  {
    _mat.resize(_n);
    for (unsigned int i = 0; i < _mat.size(); i++) {
      _mat[i].resize(_m, init);
    }
  }

  template<typename T>
  Matrix<T>::Matrix(std::vector<std::vector<T> > mat)
  : _mat(mat), _n(_mat.size()), _m(_mat[0].size()), _name(" ")
  {

  }
  template<typename T>
  Matrix<T>::Matrix(std::string name, std::vector<std::vector<T> > mat)
  : _mat(mat), _n(_mat.size()), _m(_mat[0].size()), _name(name)
  {

  }
  template<typename T>
  Matrix<T>::Matrix(unsigned int n, std::vector<T> flat)
  : _n(n), _m(n), _name(" ")
  {
    if (n*n == flat.size()) {
      _mat.resize(_n);
      for (unsigned int i = 0; i < _n; i++) {
        _mat[i].resize(_m);
        for (unsigned int j = 0; j < _m; j++) {
          this->_mat[i][j] = flat[i*_m + j];
        }
      }
    }
  }
  template<typename T>
  Matrix<T>::Matrix(std::string name, unsigned int n, std::vector<T> flat)
  : _n(n), _m(n), _name(name)
  {
    if (n*n == flat.size()) {
      _mat.resize(_n);
      for (unsigned int i = 0; i < _n; i++) {
        _mat[i].resize(_m);
        for (unsigned int j = 0; j < _m; j++) {
          this->_mat[i][j] = flat[i*_m + j];
        }
      }
    }
  }
  template<typename T>
  Matrix<T>::Matrix(unsigned int n, unsigned int m, std::vector<T> flat)
  : _n(n), _m(m), _name(" ")
  {
    if (n*m == flat.size()) {
      _mat.resize(_n);
      for (unsigned int i = 0; i < _n; i++) {
        _mat[i].resize(_m);
        for (unsigned int j = 0; j < _m; j++) {
          this->_mat[i][j] = flat[i*_m + j];
        }
      }
    }
  }

  template<typename T>
  Matrix<T>::Matrix(std::string name, unsigned int n, unsigned int m, std::vector<T> flat)
  : _n(n), _m(m), _name(name)
  {
    if (n*m == flat.size()) {
      _mat.resize(_n);
      for (unsigned int i = 0; i < _n; i++) {
        _mat[i].resize(_m);
        for (unsigned int j = 0; j < _m; j++) {
          this->_mat[i][j] = flat[i*_m + j];
        }
      }
    }
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

  //  Setters
  template<typename T>
  void Matrix<T>::set_name(std::string name)
  {
    _name = name;
  }

  //  Operator overloads
  template<typename T>
  Matrix<T>& Matrix<T>::operator=(const Matrix<T>& m)
  {
    if (&m == this)
      return *this;

    _n = m.get_rows();
    _m = m.get_cols();
    _mat.resize(_n);
    for (unsigned int i = 0; i < _n; i++) {
      _mat[i].resize(_m);
      for (unsigned int j = 0; j < _m; j++) {
        _mat[i][j] = m(i,j);
      }
    }
    return *this;
  }
  template<typename T>
  bool Matrix<T>::operator==(const Matrix<T>& m)
  {
    if (_n != m.get_rows() || _m != m.get_cols())
      return false;
    for (unsigned int i = 0; i < _n; i++) {
      for (unsigned int j = 0; j < _m; j++) {
        if (m(i,j) != this->_mat[i][j])
          return false;
      }
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
    for (unsigned int i = 0; i < _n; i++) {
      for (unsigned int j = 0; j < _m; j++) {
        l(i,j) = this->_mat[i][j] + m(i,j);
      }
    }
    return l;
  }
  template<typename T>
  Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& m)
  {
    if(_n != m.get_rows() || _m != m.get_cols())
      return *this;
    for (unsigned int i = 0; i < _n; i++) {
      for (unsigned int j = 0; j < _m; j++) {
        this->_mat[i][j] += m(i,j);
      }
    }
    return *this;
  }
  template<typename T>
  Matrix<T> Matrix<T>::operator-(const Matrix<T>& m)
  {
    if(_n != m.get_rows() || _m != m.get_cols())
      return *this;
    Matrix<T> l(_n, _m, 0.0);
    for (unsigned int i = 0; i < _n; i++) {
      for (unsigned int j = 0; j < _m; j++) {
        l(i,j) = this->_mat[i][j] - m(i,j);
      }
    }
    return l;
  }
  template<typename T>
  Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& m)
  {
    if(_n != m.get_rows() || _m != m.get_cols())
      return *this;
    for (unsigned int i = 0; i < _n; i++) {
      for (unsigned int j = 0; j < _m; j++) {
        this->_mat[i][j] -= m(i,j);
      }
    }
    return *this;
  }
  template<typename T>
  Matrix<T> Matrix<T>::operator*(const Matrix<T>& m)
  {
    if(_m != m.get_rows())
      return *this;
    Matrix<T> l(_n,m.get_cols(),0.0);
    for (unsigned int i = 0; i < _n; i++) {
      for (unsigned int j = 0; j < _m; j++) {
        for (unsigned int k = 0; k < _n; k++) {
          l(i,j) += this->_mat[i][k] * m(k,j);
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
    for (unsigned int i = 0; i < _n; i++) {
      for (unsigned int j = 0; j < _m; j++) {
        l(i,j) = this->_mat[i][j] + s;
      }
    }
    return l;
  }
  template<typename T>
  Matrix<T> Matrix<T>::operator-(const T& s)
  {
    Matrix<T> l(_n,_m,0.0);
    for (unsigned int i = 0; i < _n; i++) {
      for (unsigned int j = 0; j < _m; j++) {
        l(i,j) = this->_mat[i][j] - s;
      }
    }
    return l;
  }
  template<typename T>
  Matrix<T> Matrix<T>::operator*(const T& s)
  {
    Matrix<T> l(_n,_m,0.0);
    for (unsigned int i = 0; i < _n; i++) {
      for (unsigned int j = 0; j < _m; j++) {
        l(i,j) = this->_mat[i][j] * s;
      }
    }
    return l;
  }
  template<typename T>
  Matrix<T> Matrix<T>::operator/(const T& s)
  {
    Matrix<T> l(_n,_m,0.0);
    if (s == 0)
      return l;
    for (unsigned int i = 0; i < _n; i++) {
      for (unsigned int j = 0; j < _m; j++) {
        l(i,j) = this->_mat[i][j] + s;
      }
    }
    return l;
  }
  template<typename T>
  Matrix<T>& Matrix<T>::operator+=(const T& s)
  {
    for (unsigned int i = 0; i < _n; i++) {
      for (unsigned int j = 0; j < _m; j++) {
        this->_mat[i][j] += s;
      }
    }
    return *this;
  }
  template<typename T>
  Matrix<T>& Matrix<T>::operator-=(const T& s)
  {
    for (unsigned int i = 0; i < _n; i++) {
      for (unsigned int j = 0; j < _m; j++) {
        this->_mat[i][j] -= s;
      }
    }
    return *this;
  }
  template<typename T>
  Matrix<T>& Matrix<T>::operator*=(const T& s)
  {
    for (unsigned int i = 0; i < _n; i++) {
      for (unsigned int j = 0; j < _m; j++) {
        this->_mat[i][j] *= s;
      }
    }
    return *this;
  }
  template<typename T>
  Matrix<T>& Matrix<T>::operator/=(const T& s)
  {
    if (s == 0)
      return *this;
    for (unsigned int i = 0; i < _n; i++) {
      for (unsigned int j = 0; j < _m; j++) {
        this->_mat[i][j] += s;
      }
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
        temp += this->_mat[i][j] * v(j);
      }
      v2(i) = temp;
    }
    return v2;
  }
  //  Access Matrix<T>::operators
  template<typename T>
  T& Matrix<T>::operator()(const unsigned int& i, const unsigned int& j)
  {
    return this->_mat[i][j];
  }
  template<typename T>
  const T& Matrix<T>::operator()(const unsigned int& i, const unsigned int& j) const
  {
    return this->_mat[i][j];
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
        std::cout << this->_mat[i][j] << " ";
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
