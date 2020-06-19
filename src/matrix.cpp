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
    _array = matrix._array;
    _n = matrix.getNumRows();
    _m = matrix.getNumCols();
    _name = matrix.getName();
  }

  template<typename T>
  Matrix<T>::Matrix(uint64_t n) : _n(n), _m(n), _name(" ")
  {
    _array.resize(_n*_n,0.0);
  }

  template<typename T>
  Matrix<T>::Matrix(std::string name, uint64_t n)
  : _n(n), _m(n), _name(name)
  {
    _array.resize(_n*_n,0.0);
  }

  template<typename T>
  Matrix<T>::Matrix(uint64_t n, uint64_t m) : _n(n), _m(m), _name(" ")
  {
    _array.resize(_n*_m,0.0);
  }

  template<typename T>
  Matrix<T>::Matrix(std::string name, uint64_t n, uint64_t m)
  : _n(n), _m(m), _name(name)
  {
    _array.resize(_n*_m,0.0);
  }

  template<typename T>
  Matrix<T>::Matrix(uint64_t n, uint64_t m, const T& init)
  : _n(n), _m(m), _name(" ")
  {
    _array.resize(_n*_m, init);
  }

  template<typename T>
  Matrix<T>::Matrix(std::string name, uint64_t n,
    uint64_t m, const T& init)
  : _n(n), _m(m), _name(name)
  {
    _array.resize(_n*_m, init);
  }

  template<typename T>
  Matrix<T>::Matrix(uint64_t n, std::vector<T> flat)
  : _n(n), _m(n), _name(" "), _array(flat)
  {

  }
  template<typename T>
  Matrix<T>::Matrix(std::string name, uint64_t n, std::vector<T> flat)
  : _n(n), _m(n), _name(name), _array(flat)
  {

  }
  template<typename T>
  Matrix<T>::Matrix(uint64_t n, uint64_t m, std::vector<T> flat)
  : _n(n), _m(m), _name(" "), _array(flat)
  {

  }

  template<typename T>
  Matrix<T>::Matrix(std::string name, uint64_t n, uint64_t m, std::vector<T> flat)
  : _n(n), _m(m), _name(name), _array(flat)
  {
  }

  template<typename T>
  Matrix<T>::Matrix(std::vector<std::vector<T> > array)
  : _n(array.size()), _m(array[0].size()), _name(" ")
  {
    std::vector<T> flat;
    for (uint64_t i = 0; i < _n; i++)
    {
      flat.insert(end(flat),begin(array[i]),end(array[i]));
    }
    _array = flat;
  }

  template<typename T>
  Matrix<T>::Matrix(std::string name, std::vector<std::vector<T> > array)
  : _n(array.size()), _m(array[0].size()), _name(name)
  {
    std::vector<T> flat;
    for (uint64_t i = 0; i < _n; i++)
    {
      flat.insert(end(flat),begin(array[i]),end(array[i]));
    }
    _array = flat;
  }

  template<typename T>
  uint64_t Matrix<T>::getNumRows() const
  {
    return _n;
  }

  template<typename T>
  uint64_t Matrix<T>::getNumCols() const
  {
    return _m;
  }

  template<typename T>
  std::string Matrix<T>::getName() const
  {
    return _name;
  }

  template<typename T>
  std::vector<T> Matrix<T>::getArray() const
  {
    return _array;
  }

  template<typename T>
  float* Matrix<T>::data()
  {
    float *copy = new float[_n*_m];
    for (uint64_t i = 0; i < _n*_m; i++)
    {
      copy[i] = _array[i];
    }
    return copy;
  }

  template<typename T>
  std::vector<T> Matrix<T>::getRow(uint64_t i)
  {
    std::vector<T> row(_m,0.0);
    for (uint64_t j = 0; j < _m; j++) {
      row[j] = _array[i*_m + j];
    }
    return row;
  }

  template<typename T>
  std::vector<T> Matrix<T>::getCol(uint64_t i)
  {
    std::vector<T> col(_n,0.0);
    for (uint64_t j = 0; j < _n; j++) {
      col[j] = _array[j*_m + i];
    }
    return col;
  }

  //  Setters
  template<typename T>
  void Matrix<T>::setName(std::string name)
  {
    _name = name;
  }
  template<typename T>
  void Matrix<T>::setRow(uint64_t i, std::vector<T> row)
  {
    for (uint64_t j = 0; j < _m; j++) {
      _array[i*_m + j] = row[j];
    }
  }
  template<typename T>
  void Matrix<T>::setCol(uint64_t i, std::vector<T> col)
  {
    for (uint64_t j = 0; j < _n; j++) {
      _array[j*_m + i] = col[j];
    }
  }
  template<typename T>
  void Matrix<T>::setArray(uint64_t m, std::vector<T> mat)
  {
    _n = mat.size()/m;
    _m = m;
    _array = mat;
  }
  template<typename T>
  void Matrix<T>::setArray(std::vector<std::vector<T> > mat)
  {
    _n = mat.size();
    _m = mat[0].size();
    _array.resize(_n*_m);
    for (uint64_t i = 0; i < _n; i++)
    {
      for (uint64_t j = 0; j < _m; j++)
      {
        _array[i*_m + j] = mat[i][j];
      }
    }
  }

  //  Operator overloads
  template<typename T>
  Matrix<T>& Matrix<T>::operator=(const Matrix<T>& matrix)
  {
    if (&matrix == this)
      return *this;

    _n = matrix.getNumRows();
    _m = matrix.getNumCols();
    _name = matrix.getName();
    _array.resize(_n*_m);
    for (uint64_t i = 0; i < _n*_m; i++) {
        _array[i] = matrix(i);
    }
    return *this;
  }
  template<typename T>
  bool Matrix<T>::operator==(const Matrix<T>& matrix) const
  {
    if (_n != matrix.getNumRows() || _m != matrix.getNumCols())
      return false;
    for (uint64_t i = 0; i < _n*_m; i++) {
        if (matrix(i) != _array[i])
          return false;
    }
    return true;
  }
  template<typename T>
  bool Matrix<T>::operator!=(const Matrix<T>& matrix) const
  {
    if (_n != matrix.getNumRows() || _m != matrix.getNumCols())
      return true;
    for (uint64_t i = 0; i < _n*_m; i++) {
        if (matrix(i) != _array[i])
          return true;
    }
    return false;
  }
  template<typename T>
  Matrix<T> Matrix<T>::operator-() const
  {
    std::vector<T> mat(_n*_m);
    for (uint64_t i = 0; i < _n*_m; i++)
    {
      mat[i] = -1*_array[i];
    }
    return Matrix<T>(_name,_n,_m,mat);
  }
  //  Matrix algebra
  template<typename T>
  Matrix<T> Matrix<T>::operator+(const Matrix<T>& matrix) const
  {
    if(_n != matrix.getNumRows() || _m != matrix.getNumCols())
    {
      std::cout << "Matrices incompatible!" << std::endl;
      return *this;
    }
    std::string name = "(" + _name + " + " + matrix.getName() + ")";
    Matrix<T> l(name, _n, _m, 0.0);
    for (uint64_t i = 0; i < _n*_m; i++) {
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
    for (uint64_t i = 0; i < _n*_m; i++) {
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
    Matrix<T> l(name,_n, _m, 0.0);
    for (uint64_t i = 0; i < _n*_m; i++) {
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
    for (uint64_t i = 0; i < _n*_m; i++) {
        _array[i] -= matrix(i);
    }
    return *this;
  }
  template<typename T>
  Matrix<T> Matrix<T>::operator*(const Matrix<T>& matrix) const
  {
    if(_m != matrix.getNumRows())
    {
      std::cout << "Matrices incompatible!" << std::endl;
      return *this;
    }
    std::string name = "(" + _name + " * " + matrix.getName() + ")";
    std::vector<T> l(_n*matrix.getNumCols(),0.0);
    //  CBLAS function for matrix multiplication, A*B = C.
    /*  clbas_dgemm(Order  - either CblasRowMajor or CblasColumnMajor
                    TransA - either transpose or not for Matrix A
                    TransB - either transpose or not for Matrix B
                    M      - number of rows in A and C
                    N      - number of columns in B and C
                    K      - number of columns in A, number of rows in B
                    alpha  - scaling factor for A and B
                    A      - Matrix A, i.e. ref. to the beginning of the vector
                    Ida    - number of columns of A; _m
                    B      - Matrix B, ref. to the beginning of B
                    Idb    - number of columns of B; m.getNumCols()
                    beta   - scaling factor for C
                    C      - Matrix C, ref. to beginning of C
                    Idc    - number of columns of C; _n
        Since we are using std::vector, we must pass in a reference to
        the pointer given by .begin().
    */
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                _n, matrix.getNumCols(), _m, 1.0,
                & *_array.begin(), _m,
                & *matrix.getArray().begin(), _n,
                0.0, & *l.begin(), _n);
    return Matrix<T>(name,_n,matrix.getNumCols(),l);
  }
  template<typename T>
  Matrix<T> Matrix<T>::brute_mul(const Matrix<T>& matrix) const
  {
    if(_m != matrix.getNumRows())
    {
      std::cout << "Matrices incompatible!" << std::endl;
      return *this;
    }
    std::string name = "(" + _name + " * " + matrix.getName() + ")";
    Matrix<T> l(name,_n,matrix.getNumCols(),0.0);
    for (uint64_t i = 0; i < _n; i++) {
      for (uint64_t j = 0; j < _m; j++) {
        for (uint64_t k = 0; k < _n; k++) {
          l(i,j) += this->_array[i*_m + k] * matrix(k,j);
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
    for (uint64_t i = 0; i < _n*_m; i++) {
        l(i) = _array[i] + s;
    }
    return l;
  }
  template<typename T>
  Matrix<T> Matrix<T>::operator-(const T& s) const
  {
    std::string name = "(" + _name + " - " + std::to_string(s) + "I)";
    Matrix<T> l(name,_n,_m,0.0);
    for (uint64_t i = 0; i < _n*_m; i++) {
        l(i) = _array[i] - s;
    }
    return l;
  }
  template<typename T>
  Matrix<T> Matrix<T>::operator*(const T& s) const
  {
    std::string name = "(" + _name + " * " + std::to_string(s) + ")";
    Matrix<T> l(name,_n,_m,0.0);
    for (uint64_t i = 0; i < _n*_m; i++) {
        l(i) = _array[i] * s;
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
    for (uint64_t i = 0; i < _n*_m; i++) {
        l(i) = _array[i] + s;
    }
    return l;
  }

  template<typename T>
  Matrix<T>& Matrix<T>::operator+=(const T& s)
  {
    std::string name = "(" + _name + " + " + std::to_string(s) + "I)";
    setName(name);
    for (uint64_t i = 0; i < _n*_m; i++) {
        _array[i] += s;
    }
    return *this;
  }
  template<typename T>
  Matrix<T>& Matrix<T>::operator-=(const T& s)
  {
    std::string name = "(" + _name + " - " + std::to_string(s) + "I)";
    setName(name);
    for (uint64_t i = 0; i < _n*_m; i++) {
        _array[i] -= s;
    }
    return *this;
  }
  template<typename T>
  Matrix<T>& Matrix<T>::operator*=(const T& s)
  {
    std::string name = "(" + _name + " * " + std::to_string(s) + ")";
    setName(name);
    for (uint64_t i = 0; i < _n*_m; i++) {
        _array[i] *= s;
    }
    return *this;
  }
  template<typename T>
  Matrix<T>& Matrix<T>::operator/=(const T& s)
  {
    if (s == 0)
      return *this;
    std::string name = "(" + _name + " / " + std::to_string(s) + ")";
    setName(name);
    for (uint64_t i = 0; i < _n*_m; i++) {
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
    if (_m != v.get_dim())
      return v2;
    for (uint64_t i = 0; i < _n; i++) {
      T temp = 0.0;
      for (uint64_t j = 0; j < _m; j++) {
        temp += this->_array[i*_m + j] * v(j);
      }
      v2(i) = temp;
    }
    return v2;
  }
  //  Access Matrix<T>::operators
  template<typename T>
  T& Matrix<T>::operator()(const uint64_t& i, const uint64_t& j)
  {
    return this->_array[i*_m + j];
  }
  template<typename T>
  const T& Matrix<T>::operator()(const uint64_t& i, const uint64_t& j) const
  {
    return this->_array[i*_m + j];
  }
  template<typename T>
  T& Matrix<T>::operator()(const uint64_t& i)
  {
    return this->_array[i];
  }
  template<typename T>
  const T& Matrix<T>::operator()(const uint64_t& i) const
  {
    return this->_array[i];
  }
  template<typename T>
  std::ostream& operator<<(std::ostream& os, const Matrix<T>& matrix)
  {
    os << "(" << matrix.getNumRows() << " x " << matrix.getNumCols() << ") Matrix";
    if (matrix.getName() != " ")
      os << ": '" << matrix.getName() << "'";
    os << "\n[ ";
    for (uint64_t i = 0; i < matrix.getNumRows(); i++) {
      os << "[ ";
      for (uint64_t j = 0; j < matrix.getNumCols(); j++) {
        os << matrix(i,j) << " ";
      }
      os << "]";
      if (i < matrix.getNumRows()-1)
        os << "\n  ";
    }
    os << " ]" << std::endl;
    return os;
  }

  template<typename T>
  Matrix<T> Matrix<T>::permutationMatrix(int& n, int* pivot)
  {
    //  generate a permutation matrix from a set of pivot indices
    std::vector<T> swaps(n*n, 0);
    std::vector<uint64_t> p(n,0);
    for (int i = 0; i < n; i++)
    {
      p[i] = i;
    }
    for (int i = 0; i < n; i++)
    {
      int temp;
      temp = p[pivot[i]-1];
      p[pivot[i]-1] = p[i];
      p[i] = temp;
    }
    for(uint64_t i = 0; i < n; i++)
    {
      swaps[p[i]*n + i] = 1;
    }
    std::string name = "(" + std::to_string(n) = "x" + std::to_string(n) + ") perm";
    return Matrix<T>(name,n,n,swaps);
  }

  template<typename T>
  Matrix<T> Matrix<T>::getL(const Matrix<T>& perm)
  {
    return *this;
  }

  template<typename T>
  Matrix<T> Matrix<T>::getU(const Matrix<T>& perm)
  {
    return *this;
  }

  template<typename T>
  Matrix<T> Matrix<T>::inverse()
  {
    std::vector<T> _array_copy = _array;
    //  pivot array with indices 1 <= i <= min(n,m)
    int *ipiv = new int[_n+1];
    //  workspaces for inversion
    int lWork = _n*_m;
    double *work = new double[lWork];
    int n = _n;
    int m = _m;
    int info;
    //  first construct an LU factorization to generate the pivot indices
    //  in ipiv.
    dgetrf_(&n,&m,& *_array_copy.begin(),&n,ipiv,&info);
    dgetri_(&n,& *_array_copy.begin(),&n,ipiv,work,&lWork,&info);
    std::string name = "(" + _name + ")^-1";
    Matrix<T> inv(name,_m,_n,_array_copy);

    delete ipiv;
    delete work;
    return inv;
  }
  template<typename T>
  std::tuple<Matrix<T>,Matrix<T>,Matrix<T>> Matrix<T>::LU()
  {
    std::vector<T> _array_copy = _array;
        //  pivot array with indices 1 <= i <= min(n,m)
    int *ipiv = new int[_n+1];
    int n = _n;
    int m = _m;
    int info;
    //  first construct an LU factorization to generate the pivot indices
    //  in ipiv.
    dgetrf_(&n,&m,&*_array_copy.begin(),&n,ipiv,&info);
    Matrix<T> p = permutationMatrix(n,ipiv);
    Matrix<T> l = getL(p);
    Matrix<T> u = getU(p);
    std::tuple<Matrix<T>,Matrix<T>,Matrix<T>> result(p,l,u);
    return result;

  }
  template<typename T>
  void Matrix<T>::QR()
  {

  }
  template<typename T>
  void Matrix<T>::SVD()
  {

  }
  template<typename T>
  std::vector<T> Matrix<T>::singularValues(int info)
  {

  }

  //  Various methods
  template<typename T>
  void Matrix<T>::print()
  {
    std::cout << "(" << _n << " x " << _m << ") Matrix";
    if (_name != " ")
      std::cout << ": '" << _name << "'";
    std:: cout << "\n[ ";
    for (uint64_t i = 0; i < _n; i++) {
      for (uint64_t j = 0; j < _m; j++) {
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
    std::string sum = "dim: (" + std::to_string(_n)
                    + "x" + std::to_string(_m)
                    + "), type: "
                    + type_name<decltype(_array[0])>();
    if (_name != " ")
    {
      sum +=  ", name: '" + _name + "'";
    }
    sum += "\n[ ";
    if (_n < 10)
    {
      for (uint64_t i = 0; i < _n; i++)
      {
        if (_m < 10)
        {
          if (_array[i*_m] >= 0.0)
            sum += " ";
          for (uint64_t j = 0; j < _m; j++)
          {
            //s << std::fixed << std::setprecision(3) << this->_array[i*_m + j];
            sum += scientific_not(this->_array[i*_m + j],3);
            if (j < _m-1)
            {
              if (_array[i*_m + j+1] >= 0.0)
                sum += "   ";
              else
                sum += "  ";
            }
          }
        }
        else
        {
          if (_array[i*_m] >= 0.0)
            sum += " ";
          //s << std::fixed << std::setprecision(3) << this->_array[i*_m + 0];
          sum += scientific_not(this->_array[i*_m + 0],3);
          if (_array[i*_m + 1] >= 0.0)
            sum += "   ";
          else
            sum += "  ";
          //s << std::fixed << std::setprecision(3) << this->_array[i*_m + 1];
          sum += scientific_not(this->_array[i*_m + 1],3);
          if (_array[i*_m + 2] >= 0.0)
            sum += "   ";
          else
            sum += "  ";
          //s << std::fixed << std::setprecision(3) << this->_array[i*_m + 2];
          sum += scientific_not(this->_array[i*_m + 2],3);
          sum += "   ";
          sum += "...   ...   ...   ";
          if (_array[i*_m + _m-3] >= 0.0)
            sum += " ";
          //s << std::fixed << std::setprecision(3) << this->_array[i*_m + _m-3];
          sum += scientific_not(this->_array[i*_m + _m-3],3);
          if (_array[i*_m + _m-2] >= 0.0)
            sum += "   ";
          else
            sum += "  ";
          //s << std::fixed << std::setprecision(3) << this->_array[i*_m + _m-2];
          sum += scientific_not(this->_array[i*_m + _m-2],3);
          if (_array[i*_m + _m-1] >= 0.0)
            sum += "   ";
          else
            sum += "  ";
          //s << std::fixed << std::setprecision(3) << this->_array[i*_m + _m-1];
          sum += scientific_not(this->_array[i*_m + _m-1],3);
        }
        if (i < _n-1)
        {
          sum += "\n  ";
        }
      }
    }
    else
    {
      for (uint64_t i = 0; i < 3; i++)
      {
        if (_m < 10)
        {
          if (_array[i*_m] >= 0.0)
            sum += " ";
          for (uint64_t j = 0; j < _m; j++)
          {
            //s << std::fixed << std::setprecision(3) << this->_array[i*_m + j];
            sum += scientific_not(this->_array[i*_m + j],3);
            if (j < _m-1)
            {
              if (_array[i*_m + j+1] >= 0.0)
                sum += "   ";
              else
                sum += "  ";
            }
          }
        }
        else
        {
          if (_array[i*_m] >= 0.0)
            sum += " ";
          //s << std::fixed << std::setprecision(3) << this->_array[i*_m + 0];
          sum += scientific_not(this->_array[i*_m + 0],3);
          if (_array[i*_m + 1] >= 0.0)
            sum += "   ";
          else
            sum += "  ";
          //s << std::fixed << std::setprecision(3) << this->_array[i*_m + 1];
          sum += scientific_not(this->_array[i*_m + 1],3);
          if (_array[i*_m + 2] >= 0.0)
            sum += "   ";
          else
            sum += "  ";
          //s << std::fixed << std::setprecision(3) << this->_array[i*_m + 2];
          sum += scientific_not(this->_array[i*_m + 2],3);
          sum += "   ";
          sum += "...   ...   ...   ";
          if (_array[i*_m + _m-3] >= 0.0)
            sum += " ";
          //s << std::fixed << std::setprecision(3) << this->_array[i*_m + _m-3];
          sum += scientific_not(this->_array[i*_m + _m-3],3);
          if (_array[i*_m + _m-2] >= 0.0)
            sum += "   ";
          else
            sum += "  ";
          //s << std::fixed << std::setprecision(3) << this->_array[i*_m + _m-2];
          sum += scientific_not(this->_array[i*_m + _m-2],3);
          if (_array[i*_m + _m-1] >= 0.0)
            sum += "   ";
          else
            sum += "  ";
          //s << std::fixed << std::setprecision(3) << this->_array[i*_m + _m-1];
          sum += scientific_not(this->_array[i*_m + _m-1],3);
        }
        if (i < _n-1)
          sum += "\n  ";
      }
      sum += "    ...";
      if (_m < 10)
      {
        for (uint64_t j = 0; j < _m-1; j++)
        {
          sum += "         ...";
        }
        sum += "\n  ";
      }
      else
      {
        sum += "         ...         ...      ...   ...   ...       ...         ...         ...\n  ";
      }
      for (uint64_t i = _n-3; i < _n; i++)
      {
        if (_m < 10)
        {
          if (_array[i*_m] >= 0.0)
            sum += " ";
          for (uint64_t j = 0; j < _m; j++)
          {
            //s << std::fixed << std::setprecision(3) << this->_array[i*_m + j];
            sum += scientific_not(this->_array[i*_m + j],3);
            if (j < _m-1)
            {
              if (_array[i*_m + j+1] >= 0.0)
                sum += "   ";
              else
                sum += "  ";
            }
          }
        }
        else
        {
          if (_array[i*_m] > 0.0)
            sum += " ";
          //s << std::fixed << std::setprecision(3) << this->_array[i*_m + 0];
          sum += scientific_not(this->_array[i*_m + 0],3);
          if (_array[i*_m + 1] >= 0.0)
            sum += "   ";
          else
            sum += "  ";
          //s << std::fixed << std::setprecision(3) << this->_array[i*_m + 1];
          sum += scientific_not(this->_array[i*_m + 1],3);
          if (_array[i*_m + 2] >= 0.0)
            sum += "   ";
          else
            sum += "  ";
          //s << std::fixed << std::setprecision(3) << this->_array[i*_m + 2];
          sum += scientific_not(this->_array[i*_m + 2],3);
          sum += "   ";
          sum += "...   ...   ...   ";
          if (_array[i*_m + _m-3] >= 0.0)
            sum += " ";
          //s << std::fixed << std::setprecision(3) << this->_array[i*_m + _m-3];
          sum += scientific_not(this->_array[i*_m + _m-3],3);
          if (_array[i*_m + _m-2] >= 0.0)
            sum += "   ";
          else
            sum += "  ";
          //s << std::fixed << std::setprecision(3) << this->_array[i*_m + _m-2];
          sum += scientific_not(this->_array[i*_m + _m-2],3);
          if (_array[i*_m + _m-1] >= 0.0)
            sum += "   ";
          else
            sum += "  ";
          //s << std::fixed << std::setprecision(3) << this->_array[i*_m + _m-1];
          sum += scientific_not(this->_array[i*_m + _m-1],3);
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
  Matrix<T> identity(uint64_t n)
  {
    std::string name = std::to_string(n) + " x " + std::to_string(n) +
                        " Identity Matrix";
    Matrix<T> matrix(name, n, n, 0.0);

    for (uint64_t i = 0; i < n; i++)
    {
      matrix(i,i) = 1.0;
    }
    return matrix;
  }
  template<typename T>
  Matrix<T> zeroes(uint64_t n)
  {
    Matrix<T> z(n, n, 0.0);
    return z;
  }
  template<typename T>
  Matrix<T> zeroes(uint64_t n, uint64_t m)
  {
    Matrix<T> z(n, m, 0.0);
    return z;
  }
  template<typename T>
  Matrix<T> ones(uint64_t n)
  {
    Matrix<T> o(n, n, 1.0);
    return o;
  }
  template<typename T>
  Matrix<T> ones(uint64_t n, uint64_t m)
  {
    Matrix<T> o(n, m, 1.0);
    return o;
  }
}
