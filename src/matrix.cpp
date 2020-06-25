//  Matrix wrapper cpp file
#include "vector.h"
#include "matrix.h"

namespace ET
{

  template<typename T>
  Matrix<T>::Matrix() : _name(" "), _m(0), _n(0) {}

  template<typename T>
  Matrix<T>::~Matrix() {}

  //  Copy constructor
  template<typename T>
  Matrix<T>::Matrix(const Matrix<T>& matrix)
  {
    _array = matrix._array;
    _m = matrix.getNumRows();
    _n = matrix.getNumCols();
    _name = matrix.getName();
  }

  template<typename T>
  Matrix<T>::Matrix(uint32_t m) : _m(m), _n(m), _name(" ")
  {
    _array.resize(_m*_n,0.0);
  }

  template<typename T>
  Matrix<T>::Matrix(std::string name, uint32_t m)
  : _m(m), _n(m), _name(name)
  {
    _array.resize(_m*_n,0.0);
  }

  template<typename T>
  Matrix<T>::Matrix(uint32_t m, uint32_t n) : _m(m), _n(n), _name(" ")
  {
    _array.resize(_m*_n,0.0);
  }

  template<typename T>
  Matrix<T>::Matrix(std::string name, uint32_t m, uint32_t n)
  : _m(m), _n(n), _name(name)
  {
    _array.resize(_m*_n,0.0);
  }

  template<typename T>
  Matrix<T>::Matrix(uint32_t m, uint32_t n, const T& init)
  : _m(m), _n(n), _name(" ")
  {
    _array.resize(_m*_n, init);
  }

  template<typename T>
  Matrix<T>::Matrix(std::string name, uint32_t m,
    uint32_t n, const T& init)
  : _m(m), _n(n), _name(name)
  {
    _array.resize(_m*_n, init);
  }

  template<typename T>
  Matrix<T>::Matrix(uint32_t m, std::vector<T> flat)
  : _m(m), _n(m), _name(" "), _array(flat)
  {

  }
  template<typename T>
  Matrix<T>::Matrix(std::string name, uint32_t m, std::vector<T> flat)
  : _m(m), _n(m), _name(name), _array(flat)
  {

  }
  template<typename T>
  Matrix<T>::Matrix(uint32_t m, uint32_t n, std::vector<T> flat)
  : _m(m), _n(n), _name(" "), _array(flat)
  {

  }

  template<typename T>
  Matrix<T>::Matrix(std::string name, uint32_t m, uint32_t n, std::vector<T> flat)
  : _m(m), _n(n), _name(name), _array(flat)
  {
  }

  template<typename T>
  Matrix<T>::Matrix(std::string name, uint32_t m, uint32_t n, T* array)
  : _m(m), _n(n), _name(name)
  {
    std::vector<T> flat(array, array + _m*_n);
    _array = flat;
  }

  template<typename T>
  Matrix<T>::Matrix(std::vector<std::vector<T> > array)
  : _m(array.size()), _n(array[0].size()), _name(" ")
  {
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
    std::vector<T> flat;
    for (uint32_t i = 0; i < _m; i++)
    {
      flat.insert(end(flat),begin(array[i]),end(array[i]));
    }
    _array = flat;
  }

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

  template<typename T>
  std::vector<T> Matrix<T>::getArray() const
  {
    return _array;
  }

  template<typename T>
  std::vector<T> Matrix<T>::accessArray()
  {
    return _array;
  }

  template<typename T>
  float* Matrix<T>::data()
  {
    float *copy = new float[_m*_n];
    for (uint32_t i = 0; i < _m*_n; i++)
    {
      copy[i] = _array[i];
    }
    return copy;
  }

  template<typename T>
  std::vector<T> Matrix<T>::getRow(uint32_t i)
  {
    std::vector<T> row(_m,0.0);
    //  check that row exists
    if (i >= _m)
    {
      std::cout << "ERROR! Index " << std::to_string(i) <<
        " exceeds matrix with dimension " + std::to_string(_m) + "!";
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
    for (uint32_t j = 0; j < _n; j++) {
      col[j] = _array[j*_n + i];
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
  void Matrix<T>::setRow(uint32_t i, std::vector<T> row)
  {
    for (uint32_t j = 0; j < _n; j++) {
      _array[i*_n + j] = row[j];
    }
  }
  template<typename T>
  void Matrix<T>::setCol(uint32_t i, std::vector<T> col)
  {
    for (uint32_t j = 0; j < _m; j++) {
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

  //  Operator overloads
  template<typename T>
  Matrix<T>& Matrix<T>::operator=(const Matrix<T>& matrix)
  {
    if (&matrix == this)
      return *this;

    _m = matrix.getNumRows();
    _n = matrix.getNumCols();
    _name = matrix.getName();
    _array.resize(_m*_n);
    for (uint32_t i = 0; i < _m*_n; i++) {
        _array[i] = matrix(i);
    }
    return *this;
  }
  template<typename T>
  bool Matrix<T>::operator==(const Matrix<T>& matrix) const
  {
    if (_n != matrix.getNumRows() || _m != matrix.getNumCols())
      return false;
    for (uint32_t i = 0; i < _m*_n; i++) {
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
    for (uint32_t i = 0; i < _m*_n; i++) {
        if (matrix(i) != _array[i])
          return true;
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
    if(_n != matrix.getNumRows() || _m != matrix.getNumCols())
    {
      std::cout << "Matrices incompatible!" << std::endl;
      return *this;
    }
    std::string name = "(" + _name + " + " + matrix.getName() + ")";
    Matrix<T> l(name, _m, _n, 0.0);
    for (uint32_t i = 0; i < _m*_n; i++) {
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
    for (uint32_t i = 0; i < _m*_n; i++) {
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
    for (uint32_t i = 0; i < _m*_n; i++) {
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
    for (uint32_t i = 0; i < _m*_n; i++) {
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
    std::string name = "(" + _name + " * " + matrix.getName() + ")";
    std::vector<T> l(_m*matrix.getNumCols(),0.0);
    //  CBLAS function for matrix multiplication, A*B = C.
    /*  clbas_dgemm(Order  - either CblasRowMajor or CblasColumnMajor
                    TransA - either transpose or not for Matrix A
                    TransB - either transpose or not for Matrix B
                    M      - number of rows in A and C
                    N      - number of columns in B and C
                    K      - number of columns in A, number of rows in B
                    alpha  - scaling factor for A and B
                    A      - Matrix A, i.e. ref. to the beginning of the vector
                    Ida    - number of columns of A; _l
                    B      - Matrix B, ref. to the beginning of B
                    Idb    - number of columns of B; m.getNumCols()
                    beta   - scaling factor for C
                    C      - Matrix C, ref. to beginning of C
                    Idc    - number of columns of C; _n
        Since we are using std::vector, we must pass in a reference to
        the pointer given by .begin().
    */
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                _m, matrix.getNumCols(), _n, 1.0,
                & *_array.begin(), _n,
                & *matrix.getArray().begin(), matrix.getNumCols(),
                0.0, & *l.begin(), matrix.getNumCols());
    return Matrix<T>(name,_m,matrix.getNumCols(),l);
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
    for (uint32_t i = 0; i < _m*_n; i++) {
        l(i) = _array[i] + s;
    }
    return l;
  }
  template<typename T>
  Matrix<T> Matrix<T>::operator-(const T& s) const
  {
    std::string name = "(" + _name + " - " + std::to_string(s) + "I)";
    Matrix<T> l(name,_m,_n,0.0);
    for (uint32_t i = 0; i < _m*_n; i++) {
        l(i) = _array[i] - s;
    }
    return l;
  }
  template<typename T>
  Matrix<T> Matrix<T>::operator*(const T& s) const
  {
    std::string name = "(" + _name + " * " + std::to_string(s) + ")";
    Matrix<T> l(name,_m,_n,0.0);
    for (uint32_t i = 0; i < _m*_n; i++) {
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
      return l;
    for (uint32_t i = 0; i < _m*_n; i++) {
        l(i) = _array[i] + s;
    }
    return l;
  }

  template<typename T>
  Matrix<T>& Matrix<T>::operator+=(const T& s)
  {
    std::string name = "(" + _name + " + " + std::to_string(s) + "I)";
    setName(name);
    for (uint32_t i = 0; i < _m*_n; i++) {
        _array[i] += s;
    }
    return *this;
  }
  template<typename T>
  Matrix<T>& Matrix<T>::operator-=(const T& s)
  {
    std::string name = "(" + _name + " - " + std::to_string(s) + "I)";
    setName(name);
    for (uint32_t i = 0; i < _m*_n; i++) {
        _array[i] -= s;
    }
    return *this;
  }
  template<typename T>
  Matrix<T>& Matrix<T>::operator*=(const T& s)
  {
    std::string name = "(" + _name + " * " + std::to_string(s) + ")";
    setName(name);
    for (uint32_t i = 0; i < _m*_n; i++) {
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
    for (uint32_t i = 0; i < _m*_n; i++) {
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
      return v2;
    for (uint32_t i = 0; i < _m; i++) {
      T temp = 0.0;
      for (uint32_t j = 0; j < _n; j++) {
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
    return this->_array[i*_n + j];
  }
  template<typename T>
  const T& Matrix<T>::operator()(const uint32_t& i, const uint32_t& j) const
  {
    return this->_array[i*_n + j];
  }
  template<typename T>
  T& Matrix<T>::operator()(const uint32_t& i)
  {
    return this->_array[i];
  }
  template<typename T>
  const T& Matrix<T>::operator()(const uint32_t& i) const
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
    for (uint32_t i = 0; i < matrix.getNumRows(); i++) {
      os << "[ ";
      for (uint32_t j = 0; j < matrix.getNumCols(); j++) {
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
    std::vector<uint32_t> p(n,0);
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
    for(uint32_t i = 0; i < n; i++)
    {
      swaps[p[i]*n + i] = 1;
    }
    std::string name = "(" + std::to_string(n) = "x" + std::to_string(n) + ") perm";
    return Matrix<T>(name,n,n,swaps);
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
  Matrix<T> Matrix<T>::inverse()
  {
    int m = _m;
    int n = _n;
    std::vector<T> _array_copy = _array;
    //  pivot array with indices 1 <= i <= min(n,m)
    int *pivot = new int[_n+1];
    //  workspaces for inversion
    int lWork = _m*_n;
    double *work = new double[_m*_n];
    int info;
    //  first construct an LU factorization to generate the pivot indices
    //  in pivot.
    dgetrf_(&m,&n,& *_array_copy.begin(),&m,pivot,&info);
    dgetri_(&m,& *_array_copy.begin(),&m,pivot,work,&lWork,&info);
    std::string name = "(" + _name + ")^-1";
    Matrix<T> inv(name,m,n,_array_copy);
    return inv;
  }
  template<typename T>
  Matrix<T> Matrix<T>::pseudoInverse()
  {
    //  first compute SVD decomposition
    std::tuple<Matrix<T>,Matrix<T>,Matrix<T>> svd = SVD();
    //  The pseudo inverse is then V * Sigma * (U)^T
    Matrix<T> U_matrix = std::get<0>(svd);
    Matrix<T> Sigma_matrix = std::get<1>(svd);
    Matrix<T> VT_matrix = std::get<2>(svd);
    U_matrix.transpose_inplace();
    VT_matrix.transpose_inplace();
    for (uint32_t i = 0; i < Sigma_matrix.getNumCols(); i++)
    {
      T value = Sigma_matrix(i,i);
      if (value > 1e-10)
        Sigma_matrix(i,i) = 1/value;
      else
        Sigma_matrix(i,i) = 0;
    }
    Sigma_matrix.transpose_inplace();
    std::string name = "(" + _name + ")^+";
    Matrix<T> result = VT_matrix * (Sigma_matrix * U_matrix);
    result.setName(name);
    return result;
  }

  template<typename T>
  std::tuple<Matrix<T>,Matrix<T>,Matrix<T>> Matrix<T>::LU()
  {
    //  pivot array with indices 1 <= i <= min(n,m)
    int *pivot = new int[_n+1];
    int info;
    int m = _m;
    int n = _n;
    //  first construct an LU factorization to generate the pivot indices
    //  in pivot.
    dgetrf_(&m,&n,&*_array.begin(),&m,pivot,&info);
    Matrix<T> P = permutationMatrix(m,pivot);
    Matrix<T> L = getL(P);
    Matrix<T> U = getU(P);
    std::tuple<Matrix<T>,Matrix<T>,Matrix<T>> result(P,L,U);
    return result;

  }
  template<typename T>
  std::tuple<Matrix<T>,Matrix<T>> Matrix<T>::QR()
  {

  }
  template<typename T>
  std::tuple<Matrix<T>,Matrix<T>,Matrix<T>> Matrix<T>::SVD()
  {
    int info;
    int m = _m;
    int n = _n;
    //  Sigma, U and V_T matrices
    //  Sigma will be a vector of singular values of size n,
    //  U is a unitary m x m matrix,
    //  VT is the transpose of an n x n unitary matrix.
    //  the ldu, lda and ldvt are the "leading dimension" (i.e. the number
    //  of rows).
    T sigma[std::min(m,n)], U[m*m], VT[n*n];
    //  workspaces for inversion
    int lWork;
    //  find the optimal workspace first
    double workOpt;
    double *work;
    lWork = -1;
    dgesvd_("A", "A", &m, &n, &*_array.begin(), &m, sigma, U,
        &m, VT, &n, &workOpt, &lWork, &info);
    lWork = (int)workOpt;
    work = (double*)malloc( lWork*sizeof(double) );
    /* Compute SVD */
    dgesvd_("A", "A", &m, &n, &*_array.begin(), &m, sigma, U,
        &m, VT, &n, work, &lWork, &info);
    /* Check for convergence */
    if( info > 0 ) {
      printf( "The algorithm computing SVD failed to converge.\n" );
      exit( 1 );
    }
    std::vector<T> sig(sigma, sigma + sizeof(sigma)/sizeof(sigma[0]));
    _singular_values = sig;
    //  construct U, VT and Sigma matrices
    Matrix<T> U_latrix("U",_m,_m,U);
    Matrix<T> VT_latrix("(V)^T",_n,_n,VT);
    Matrix<T> S_latrix("Sigma",_m,_n,0.0);
    for (uint32_t i = 0; i < _singular_values.size(); i++)
    {
      S_latrix(i,i) = _singular_values[i];
    }
    //U_latrix.transpose(true);
    //U_latrix.setName("U");
    //VT_latrix.transpose(true);
    //VT_latrix.setName("(V)^T");
    std::tuple<Matrix<T>,Matrix<T>,Matrix<T>> result(U_latrix, S_latrix, VT_latrix);
    return result;
  }

  template<typename T>
  bool Matrix<T>::isInvertible()
  {

  }

  template<typename T>
  void Matrix<T>::findSingularValues()
  {
    int info;
    int m = _m;
    int n = _n;
    //  Sigma, U and V_T matrices
    //  Sigma will be a vector of singular values of size n,
    //  U is a unitary m x m matrix,
    //  VT is the transpose of an n x n unitary matrix.
    //  the ldu, lda and ldvt are the "leading dimension" (i.e. the number
    //  of rows).
    T sigma[std::min(m,n)], U[m*m], VT[n*n];
    //  workspaces for inversion
    int lWork;
    //  find the optimal workspace first
    double workOpt;
    double *work;
    lWork = -1;
    dgesvd_("A", "A", &m, &n, &*_array.begin(), &m, sigma, U,
        &m, VT, &n, &workOpt, &lWork, &info);
    lWork = (int)workOpt;
    work = (double*)malloc( lWork*sizeof(double) );
    /* Compute SVD */
    dgesvd_("A", "A", &m, &n, &*_array.begin(), &m, sigma, U,
        &m, VT, &n, work, &lWork, &info);
    /* Check for convergence */
    if( info > 0 ) {
      printf( "The algorithm computing SVD failed to converge.\n" );
      exit( 1 );
    }
    std::vector<T> sig(sigma, sigma + sizeof(sigma)/sizeof(sigma[0]));
    _singular_values = sig;
  }

  template<typename T>
  std::vector<T> Matrix<T>::getSingularValues()
  {
    if (_singular_values.size() == 0)
      findSingularValues();
    return _singular_values;
  }

  //  Various methods
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
            //s << std::fixed << std::setprecision(3) << this->_array[i*_l + j];
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
          //s << std::fixed << std::setprecision(3) << this->_array[i*_l + 0];
          sum += scientific_not(this->_array[i*_n + 0],3);
          if (_array[i*_n + 1] >= 0.0)
            sum += "   ";
          else
            sum += "  ";
          //s << std::fixed << std::setprecision(3) << this->_array[i*_l + 1];
          sum += scientific_not(this->_array[i*_n + 1],3);
          if (_array[i*_n + 2] >= 0.0)
            sum += "   ";
          else
            sum += "  ";
          //s << std::fixed << std::setprecision(3) << this->_array[i*_l + 2];
          sum += scientific_not(this->_array[i*_n + 2],3);
          sum += "   ";
          sum += "...   ...   ...   ";
          if (_array[i*_n + _n-3] >= 0.0)
            sum += " ";
          //s << std::fixed << std::setprecision(3) << this->_array[i*_l + _l-3];
          sum += scientific_not(this->_array[i*_n + _n-3],3);
          if (_array[i*_n + _n-2] >= 0.0)
            sum += "   ";
          else
            sum += "  ";
          //s << std::fixed << std::setprecision(3) << this->_array[i*_l + _l-2];
          sum += scientific_not(this->_array[i*_n + _n-2],3);
          if (_array[i*_n + _n-1] >= 0.0)
            sum += "   ";
          else
            sum += "  ";
          //s << std::fixed << std::setprecision(3) << this->_array[i*_l + _l-1];
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
            //s << std::fixed << std::setprecision(3) << this->_array[i*_l + j];
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
          //s << std::fixed << std::setprecision(3) << this->_array[i*_l + 0];
          sum += scientific_not(this->_array[i*_n + 0],3);
          if (_array[i*_n + 1] >= 0.0)
            sum += "   ";
          else
            sum += "  ";
          //s << std::fixed << std::setprecision(3) << this->_array[i*_l + 1];
          sum += scientific_not(this->_array[i*_n + 1],3);
          if (_array[i*_n + 2] >= 0.0)
            sum += "   ";
          else
            sum += "  ";
          //s << std::fixed << std::setprecision(3) << this->_array[i*_l + 2];
          sum += scientific_not(this->_array[i*_n + 2],3);
          sum += "   ";
          sum += "...   ...   ...   ";
          if (_array[i*_n + _n-3] >= 0.0)
            sum += " ";
          //s << std::fixed << std::setprecision(3) << this->_array[i*_l + _l-3];
          sum += scientific_not(this->_array[i*_n + _n-3],3);
          if (_array[i*_n + _n-2] >= 0.0)
            sum += "   ";
          else
            sum += "  ";
          //s << std::fixed << std::setprecision(3) << this->_array[i*_l + _l-2];
          sum += scientific_not(this->_array[i*_n + _n-2],3);
          if (_array[i*_n + _n-1] >= 0.0)
            sum += "   ";
          else
            sum += "  ";
          //s << std::fixed << std::setprecision(3) << this->_array[i*_l + _l-1];
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
            //s << std::fixed << std::setprecision(3) << this->_array[i*_l + j];
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
          //s << std::fixed << std::setprecision(3) << this->_array[i*_l + 0];
          sum += scientific_not(this->_array[i*_n + 0],3);
          if (_array[i*_n + 1] >= 0.0)
            sum += "   ";
          else
            sum += "  ";
          //s << std::fixed << std::setprecision(3) << this->_array[i*_l + 1];
          sum += scientific_not(this->_array[i*_n + 1],3);
          if (_array[i*_n + 2] >= 0.0)
            sum += "   ";
          else
            sum += "  ";
          //s << std::fixed << std::setprecision(3) << this->_array[i*_l + 2];
          sum += scientific_not(this->_array[i*_n + 2],3);
          sum += "   ";
          sum += "...   ...   ...   ";
          if (_array[i*_n + _n-3] >= 0.0)
            sum += " ";
          //s << std::fixed << std::setprecision(3) << this->_array[i*_l + _l-3];
          sum += scientific_not(this->_array[i*_n + _n-3],3);
          if (_array[i*_n + _n-2] >= 0.0)
            sum += "   ";
          else
            sum += "  ";
          //s << std::fixed << std::setprecision(3) << this->_array[i*_l + _l-2];
          sum += scientific_not(this->_array[i*_n + _n-2],3);
          if (_array[i*_n + _n-1] >= 0.0)
            sum += "   ";
          else
            sum += "  ";
          //s << std::fixed << std::setprecision(3) << this->_array[i*_l + _l-1];
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

  template<typename T>
  Matrix<T> identity(uint32_t m)
  {
    std::string name = "I_{" + std::to_string(m) + "x"
                     + std::to_string(m) + "}";
    Matrix<T> matrix(name, m, m, 0.0);

    for (uint32_t i = 0; i < m; i++)
    {
      matrix(i,i) = 1.0;
    }
    return matrix;
  }
  template<typename T>
  Matrix<T> zeroes(uint32_t m)
  {
    Matrix<T> z(m, m, 0.0);
    return z;
  }
  template<typename T>
  Matrix<T> zeroes(uint32_t m, uint32_t n)
  {
    Matrix<T> z(m, n, 0.0);
    return z;
  }
  template<typename T>
  Matrix<T> ones(uint32_t m)
  {
    Matrix<T> o(m, m, 1.0);
    return o;
  }
  template<typename T>
  Matrix<T> ones(uint32_t m, uint32_t n)
  {
    Matrix<T> o(m, n, 1.0);
    return o;
  }

  //  Level 2 BLAS double
  Vector<double> DGEMV(double& alpha, Matrix<double>& A,
      Vector<double>& x)
  {
    std::vector<double> y(x.getDim());
    cblas_dgemv(CblasRowMajor,CblasNoTrans,A.getNumRows(),A.getNumCols(),
      alpha,A.getArray().data(),A.getNumRows(),
      x.getVec().data(),1,1,y.data(),1);
    std::string name;
    if (alpha != 0)
    {
      name += std::to_string(alpha) + " * ";
    }
    name += A.getName() + " * " + x.getName();
    Vector<double> result(name,y);
    return result;
  }

  Vector<double> DGEMV(double& alpha, Matrix<double>& A,
      Vector<double>& x, double& beta, Vector<double>& y)
  {
    std::vector<double> y2 = y.getVec();
    cblas_dgemv(CblasRowMajor,CblasNoTrans,A.getNumRows(),A.getNumCols(),
      alpha,A.getArray().data(),A.getNumRows(),
      x.getVec().data(),1,beta,y2.data(),1);
    Vector<double> result(y2);
    return result;
  }

}
