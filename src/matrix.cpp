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
  Matrix<T>::Matrix(const Matrix<T>& matrix)
  {
    _mat = matrix._mat;
    _n = matrix.get_rows();
    _m = matrix.get_cols();
    _name = matrix.get_name();
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
  Matrix<T>::Matrix(std::vector<std::vector<T> > array)
  : _n(array.size()), _m(array[0].size()), _name(" ")
  {
    std::vector<T> flat;
    for (unsigned int i = 0; i < _n; i++)
    {
      flat.insert(end(flat),begin(array[i]),end(array[i]));
    }
    _mat = flat;
  }

  template<typename T>
  Matrix<T>::Matrix(std::string name, std::vector<std::vector<T> > array)
  : _n(array.size()), _m(array[0].size()), _name(name)
  {
    std::vector<T> flat;
    for (unsigned int i = 0; i < _n; i++)
    {
      flat.insert(end(flat),begin(array[i]),end(array[i]));
    }
    _mat = flat;
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
  template<typename T>
  void Matrix<T>::set_mat(unsigned int m, std::vector<T> mat)
  {
    _n = mat.size()/m;
    _m = m;
    _mat = mat;
  }
  template<typename T>
  void Matrix<T>::set_mat(std::vector<std::vector<T> > mat)
  {
    _n = mat.size();
    _m = mat[0].size();
    _mat.resize(_n*_m);
    for (unsigned int i = 0; i < _n; i++)
    {
      for (unsigned int j = 0; j < _m; j++)
      {
        _mat[i*_m + j] = mat[i][j];
      }
    }
  }

  //  Operator overloads
  template<typename T>
  Matrix<T>& Matrix<T>::operator=(const Matrix<T>& matrix)
  {
    if (&matrix == this)
      return *this;

    _n = matrix.get_rows();
    _m = matrix.get_cols();
    _name = matrix.get_name();
    _mat.resize(_n*_m);
    for (unsigned int i = 0; i < _n*_m; i++) {
        _mat[i] = matrix(i);
    }
    return *this;
  }
  template<typename T>
  bool Matrix<T>::operator==(const Matrix<T>& matrix) const
  {
    if (_n != matrix.get_rows() || _m != matrix.get_cols())
      return false;
    for (unsigned int i = 0; i < _n*_m; i++) {
        if (matrix(i) != _mat[i])
          return false;
    }
    return true;
  }
  template<typename T>
  bool Matrix<T>::operator!=(const Matrix<T>& matrix) const
  {
    if (_n != matrix.get_rows() || _m != matrix.get_cols())
      return true;
    for (unsigned int i = 0; i < _n*_m; i++) {
        if (matrix(i) != _mat[i])
          return true;
    }
    return false;
  }
  template<typename T>
  Matrix<T> Matrix<T>::operator-() const
  {
    std::vector<T> mat(_n*_m);
    for (unsigned int i = 0; i < _n*_m; i++)
    {
      mat[i] = -1*_mat[i];
    }
    return Matrix<T>(_name,_n,_m,mat);
  }
  //  Matrix algebra
  template<typename T>
  Matrix<T> Matrix<T>::operator+(const Matrix<T>& matrix) const
  {
    if(_n != matrix.get_rows() || _m != matrix.get_cols())
    {
      std::cout << "Matrices incompatible!" << std::endl;
      return *this;
    }
    std::string name = "(" + _name + " + " + matrix.get_name() + ")";
    Matrix<T> l(name, _n, _m, 0.0);
    for (unsigned int i = 0; i < _n*_m; i++) {
        l(i) = _mat[i] + matrix(i);
    }
    return l;
  }
  template<typename T>
  Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& matrix)
  {
    if(_n != matrix.get_rows() || _m != matrix.get_cols())
    {
      std::cout << "Matrices incompatible!" << std::endl;
      return *this;
    }
    std::string name = "(" + _name + " + " + matrix.get_name() + ")";
    set_name(name);
    for (unsigned int i = 0; i < _n*_m; i++) {
        _mat[i] += matrix(i);
    }
    return *this;
  }
  template<typename T>
  Matrix<T> Matrix<T>::operator-(const Matrix<T>& matrix) const
  {
    if(_n != matrix.get_rows() || _m != matrix.get_cols())
    {
      std::cout << "Matrices incompatible!" << std::endl;
      return *this;
    }
    std::string name = "(" + _name + " - " + matrix.get_name() + ")";
    Matrix<T> l(name,_n, _m, 0.0);
    for (unsigned int i = 0; i < _n*_m; i++) {
        l(i) = _mat[i] - matrix(i);
    }
    return l;
  }
  template<typename T>
  Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& matrix)
  {
    if(_n != matrix.get_rows() || _m != matrix.get_cols())
    {
      std::cout << "Matrices incompatible!" << std::endl;
      return *this;
    }
    std::string name = "(" + _name + " - " + matrix.get_name() + ")";
    set_name(name);
    for (unsigned int i = 0; i < _n*_m; i++) {
        _mat[i] -= matrix(i);
    }
    return *this;
  }
  template<typename T>
  Matrix<T> Matrix<T>::operator*(const Matrix<T>& matrix) const
  {
    if(_m != matrix.get_rows())
    {
      std::cout << "Matrices incompatible!" << std::endl;
      return *this;
    }
    std::string name = "(" + _name + " * " + matrix.get_name() + ")";
    std::vector<T> l(_n*matrix.get_cols(),0.0);
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
                _n, matrix.get_cols(), _m, 1.0,
                & *_mat.begin(), _n,
                & *matrix.get_mat().begin(), _m,
                0.0, & *l.begin(), _n);
    return Matrix<T>(name,_n,matrix.get_cols(),l);
  }
  template<typename T>
  Matrix<T> Matrix<T>::brute_mul(const Matrix<T>& matrix) const
  {
    if(_m != matrix.get_rows())
    {
      std::cout << "Matrices incompatible!" << std::endl;
      return *this;
    }
    std::string name = "(" + _name + " * " + matrix.get_name() + ")";
    Matrix<T> l(name,_n,matrix.get_cols(),0.0);
    for (unsigned int i = 0; i < _n; i++) {
      for (unsigned int j = 0; j < _m; j++) {
        for (unsigned int k = 0; k < _n; k++) {
          l(i,j) += this->_mat[i*_m + k] * matrix(k,j);
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
    Matrix<T> l(name,_n,_m,0.0);
    for (unsigned int i = 0; i < _n*_m; i++) {
        l(i) = _mat[i] + s;
    }
    return l;
  }
  template<typename T>
  Matrix<T> Matrix<T>::operator-(const T& s) const
  {
    std::string name = "(" + _name + " - " + std::to_string(s) + "I)";
    Matrix<T> l(name,_n,_m,0.0);
    for (unsigned int i = 0; i < _n*_m; i++) {
        l(i) = _mat[i] - s;
    }
    return l;
  }
  template<typename T>
  Matrix<T> Matrix<T>::operator*(const T& s) const
  {
    std::string name = "(" + _name + " * " + std::to_string(s) + ")";
    Matrix<T> l(name,_n,_m,0.0);
    for (unsigned int i = 0; i < _n*_m; i++) {
        l(i) = _mat[i] * s;
    }
    return l;
  }
  template<typename T>
  Matrix<T> Matrix<T>::operator/(const T& s) const
  {
    std::string name = "(" + _name + " / " + std::to_string(s) + ")";
    Matrix<T> l(name,_n,_m,0.0);
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
    std::string name = "(" + _name + " + " + std::to_string(s) + "I)";
    set_name(name);
    for (unsigned int i = 0; i < _n*_m; i++) {
        _mat[i] += s;
    }
    return *this;
  }
  template<typename T>
  Matrix<T>& Matrix<T>::operator-=(const T& s)
  {
    std::string name = "(" + _name + " - " + std::to_string(s) + "I)";
    set_name(name);
    for (unsigned int i = 0; i < _n*_m; i++) {
        _mat[i] -= s;
    }
    return *this;
  }
  template<typename T>
  Matrix<T>& Matrix<T>::operator*=(const T& s)
  {
    std::string name = "(" + _name + " * " + std::to_string(s) + ")";
    set_name(name);
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
    std::string name = "(" + _name + " / " + std::to_string(s) + ")";
    set_name(name);
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
  std::ostream& operator<<(std::ostream& os, const Matrix<T>& matrix)
  {
    os << "(" << matrix.get_rows() << " x " << matrix.get_cols() << ") Matrix";
    if (matrix.get_name() != " ")
      os << ": '" << matrix.get_name() << "'";
    os << "\n[ ";
    for (unsigned int i = 0; i < matrix.get_rows(); i++) {
      os << "[ ";
      for (unsigned int j = 0; j < matrix.get_cols(); j++) {
        os << matrix(i,j) << " ";
      }
      os << "]";
      if (i < matrix.get_rows()-1)
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
      for (unsigned int j = 0; j < _m; j++) {
        std::cout << this->_mat[i*_m + j];
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
    std::string sum = "dim: (" + std::to_string(_n)
                    + "x" + std::to_string(_m)
                    + "), type: "
                    + type_name<decltype(_mat[0])>();
    if (_name != " ")
    {
      sum +=  ", name: '" + _name + "'";
    }
    sum += "\n[ ";
    if (_n < 10)
    {
      for (unsigned int i = 0; i < _n; i++)
      {
        if (_m < 10)
        {
          if (_mat[i*_m] > 0.0)
            sum += " ";
          for (unsigned int j = 0; j < _m; j++)
          {
            //s << std::fixed << std::setprecision(3) << this->_mat[i*_m + j];
            sum += scientific_not(this->_mat[i*_m + j],3);
            if (j < _m-1)
            {
              if (_mat[i*_m + j+1] > 0.0)
                sum += "   ";
              else
                sum += "  ";
            }
          }
        }
        else
        {
          if (_mat[i*_m] > 0.0)
            sum += " ";
          //s << std::fixed << std::setprecision(3) << this->_mat[i*_m + 0];
          sum += scientific_not(this->_mat[i*_m + 0],3);
          if (_mat[i*_m + 1] > 0.0)
            sum += "   ";
          else
            sum += "  ";
          //s << std::fixed << std::setprecision(3) << this->_mat[i*_m + 1];
          sum += scientific_not(this->_mat[i*_m + 1],3);
          if (_mat[i*_m + 2] > 0.0)
            sum += "   ";
          else
            sum += "  ";
          //s << std::fixed << std::setprecision(3) << this->_mat[i*_m + 2];
          sum += scientific_not(this->_mat[i*_m + 2],3);
          sum += "   ";
          sum += "...   ...   ...   ";
          if (_mat[i*_m + _m-3] > 0.0)
            sum += " ";
          //s << std::fixed << std::setprecision(3) << this->_mat[i*_m + _m-3];
          sum += scientific_not(this->_mat[i*_m + _m-3],3);
          if (_mat[i*_m + _m-2] > 0.0)
            sum += "   ";
          else
            sum += "  ";
          //s << std::fixed << std::setprecision(3) << this->_mat[i*_m + _m-2];
          sum += scientific_not(this->_mat[i*_m + _m-2],3);
          if (_mat[i*_m + _m-1] > 0.0)
            sum += "   ";
          else
            sum += "  ";
          //s << std::fixed << std::setprecision(3) << this->_mat[i*_m + _m-1];
          sum += scientific_not(this->_mat[i*_m + _m-1],3);
        }
        if (i < _n-1)
        {
          sum += "\n  ";
        }
      }
    }
    else
    {
      for (unsigned int i = 0; i < 3; i++)
      {
        if (_m < 10)
        {
          if (_mat[i*_m] > 0.0)
            sum += " ";
          for (unsigned int j = 0; j < _m; j++)
          {
            //s << std::fixed << std::setprecision(3) << this->_mat[i*_m + j];
            sum += scientific_not(this->_mat[i*_m + j],3);
            if (j < _m-1)
            {
              if (_mat[i*_m + j+1] > 0.0)
                sum += "   ";
              else
                sum += "  ";
            }
          }
        }
        else
        {
          if (_mat[i*_m] > 0.0)
            sum += " ";
          //s << std::fixed << std::setprecision(3) << this->_mat[i*_m + 0];
          sum += scientific_not(this->_mat[i*_m + 0],3);
          if (_mat[i*_m + 1] > 0.0)
            sum += "   ";
          else
            sum += "  ";
          //s << std::fixed << std::setprecision(3) << this->_mat[i*_m + 1];
          sum += scientific_not(this->_mat[i*_m + 1],3);
          if (_mat[i*_m + 2] > 0.0)
            sum += "   ";
          else
            sum += "  ";
          //s << std::fixed << std::setprecision(3) << this->_mat[i*_m + 2];
          sum += scientific_not(this->_mat[i*_m + 2],3);
          sum += "   ";
          sum += "...   ...   ...   ";
          if (_mat[i*_m + _m-3] > 0.0)
            sum += " ";
          //s << std::fixed << std::setprecision(3) << this->_mat[i*_m + _m-3];
          sum += scientific_not(this->_mat[i*_m + _m-3],3);
          if (_mat[i*_m + _m-2] > 0.0)
            sum += "   ";
          else
            sum += "  ";
          //s << std::fixed << std::setprecision(3) << this->_mat[i*_m + _m-2];
          sum += scientific_not(this->_mat[i*_m + _m-2],3);
          if (_mat[i*_m + _m-1] > 0.0)
            sum += "   ";
          else
            sum += "  ";
          //s << std::fixed << std::setprecision(3) << this->_mat[i*_m + _m-1];
          sum += scientific_not(this->_mat[i*_m + _m-1],3);
        }
        if (i < _n-1)
          sum += "\n  ";
      }
      sum += "    ...";
      if (_m < 10)
      {
        for (unsigned int j = 0; j < _m-1; j++)
        {
          sum += "         ...";
        }
        sum += "\n  ";
      }
      else
      {
        sum += "         ...         ...      ...   ...   ...       ...         ...         ...\n  ";
      }
      for (unsigned int i = _n-3; i < _n; i++)
      {
        if (_m < 10)
        {
          if (_mat[i*_m] > 0.0)
            sum += " ";
          for (unsigned int j = 0; j < _m; j++)
          {
            //s << std::fixed << std::setprecision(3) << this->_mat[i*_m + j];
            sum += scientific_not(this->_mat[i*_m + j],3);
            if (j < _m-1)
            {
              if (_mat[i*_m + j+1] > 0.0)
                sum += "   ";
              else
                sum += "  ";
            }
          }
        }
        else
        {
          if (_mat[i*_m] > 0.0)
            sum += " ";
          //s << std::fixed << std::setprecision(3) << this->_mat[i*_m + 0];
          sum += scientific_not(this->_mat[i*_m + 0],3);
          if (_mat[i*_m + 1] > 0.0)
            sum += "   ";
          else
            sum += "  ";
          //s << std::fixed << std::setprecision(3) << this->_mat[i*_m + 1];
          sum += scientific_not(this->_mat[i*_m + 1],3);
          if (_mat[i*_m + 2] > 0.0)
            sum += "   ";
          else
            sum += "  ";
          //s << std::fixed << std::setprecision(3) << this->_mat[i*_m + 2];
          sum += scientific_not(this->_mat[i*_m + 2],3);
          sum += "   ";
          sum += "...   ...   ...   ";
          if (_mat[i*_m + _m-3] > 0.0)
            sum += " ";
          //s << std::fixed << std::setprecision(3) << this->_mat[i*_m + _m-3];
          sum += scientific_not(this->_mat[i*_m + _m-3],3);
          if (_mat[i*_m + _m-2] > 0.0)
            sum += "   ";
          else
            sum += "  ";
          //s << std::fixed << std::setprecision(3) << this->_mat[i*_m + _m-2];
          sum += scientific_not(this->_mat[i*_m + _m-2],3);
          if (_mat[i*_m + _m-1] > 0.0)
            sum += "   ";
          else
            sum += "  ";
          //s << std::fixed << std::setprecision(3) << this->_mat[i*_m + _m-1];
          sum += scientific_not(this->_mat[i*_m + _m-1],3);
        }
        if (i < _n-1)
          sum += "\n  ";
      }
    }
    sum += "  ]";
    return sum;
  }

  template<typename T>
  Matrix<T> Matrix<T>::transpose()
  {

  }

  template<typename T>
  Matrix<T> identity(unsigned int n)
  {
    std::string name = std::to_string(n) + " x " + std::to_string(n) +
                        " Identity Matrix";
    Matrix<T> matrix(name, n, n, 0.0);

    for (unsigned int i = 0; i < n; i++)
    {
      matrix(i,i) = 1.0;
    }
    return matrix;
  }
  template<typename T>
  Matrix<T> zeroes(unsigned int n)
  {
    Matrix<T> z(n, n, 0.0);
    return z;
  }
  template<typename T>
  Matrix<T> zeroes(unsigned int n, unsigned int m)
  {
    Matrix<T> z(n, m, 0.0);
    return z;
  }
  template<typename T>
  Matrix<T> ones(unsigned int n)
  {
    Matrix<T> o(n, n, 1.0);
    return o;
  }
  template<typename T>
  Matrix<T> ones(unsigned int n, unsigned int m)
  {
    Matrix<T> o(n, m, 1.0);
    return o;
  }
}
