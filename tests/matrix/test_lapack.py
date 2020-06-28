import numpy as np
import matplotlib.pyplot as plt
from etraj.etraj import Matrix, Vector
from scipy.linalg import blas
import etraj.etraj as et

#  tests for LAPACK methods
print("\n----------Tests for LAPACK functions-----------")
print("\n-----------------------------------------------")

print("\n(1) DGELS")
'''
    Example from: https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_lapack_examples/dgels_ex.c.htm

    DGELS Example.
    ==============

    Program computes the least squares solution to the overdetermined linear
    system A*X = B with full rank matrix A using QR factorization,
    where A is the coefficient matrix:

     1.44  -7.84  -4.39   4.53
    -9.96  -0.28  -3.24   3.83
    -7.55   3.24   6.27  -6.64
     8.34   8.09   5.28   2.06
     7.08   2.52   0.74  -2.47
    -5.45  -5.70  -1.19   4.70

    and B is the right-hand side matrix:

     8.58   9.35
     8.26  -4.43
     8.48  -0.70
    -5.28  -0.26
     5.72  -7.36
     8.93  -2.52

    Description.
    ============

    The routine solves overdetermined or underdetermined real linear systems
    involving an m-by-n matrix A, or its transpose, using a QR or LQ
    factorization of A. It is assumed that A has full rank.

    Several right hand side vectors b and solution vectors x can be handled
    in a single call; they are stored as the columns of the m-by-nrhs right
    hand side matrix B and the n-by-nrhs solution matrix X.

    Example Program Results.
    ========================

    DGELS Example Program Results

    Solution
    -0.45   0.25
    -0.85  -0.90
    0.71   0.63
    0.13   0.14
'''
print("\n(1a) - DGELS")
print("Solve the linear system ac = b, with\n")
print(">>> a = [[ 1.44, -7.84, -4.39,  4.53],\n"
    + "         [-9.96, -0.28, -3.24,  3.83],\n"
    + "         [-7.55,  3.24,  6.27, -6.64],\n"
    + "         [ 8.34,  8.09,  5.28,  2.06],\n"
    + "         [ 7.08,  2.52,  0.74, -2.47],\n"
    + "         [-5.45, -5.70, -1.19,  4.70]]")
print(">>> b = [[ 8.58,  9.35],\n"
    + "         [ 8.26, -4.43],\n"
    + "         [ 8.48, -0.70],\n"
    + "         [-5.28, -0.26],\n"
    + "         [ 5.72, -7.36],\n"
    + "         [ 8.93, -2.52]]")
print(">>> A = Matrix('A',a)")
print(">>> B = Matrix('B',b)")
print(">>> C = et.dgels(A,B)")
print(">>> print(C)")
a = [[1.44, -7.84, -4.39,  4.53],
    [-9.96, -0.28, -3.24,  3.83],
    [-7.55,  3.24,  6.27, -6.64],
    [ 8.34,  8.09,  5.28,  2.06],
    [ 7.08,  2.52,  0.74, -2.47],
    [-5.45, -5.70, -1.19,  4.70]]
b = [[ 8.58,  9.35],
     [ 8.26, -4.43],
     [ 8.48, -0.70],
     [-5.28, -0.26],
     [ 5.72, -7.36],
     [ 8.93, -2.52]]
A = Matrix('A',a)
B = Matrix('B',b)
C = et.dgels(A,B)
print(C)

'''
    Example from: http://www.netlib.org/lapack/explore-html/d6/db6/a01172_source.html

    Description
    ===========

    In this example, we wish solve the least squares problem min_x || B - Ax ||
    for two right-hand sides using the LAPACK routine DGELS. For input we will
    use the 5-by-3 matrix

          ( 1  1  1 )
          ( 2  3  4 )
      A = ( 3  5  2 )
          ( 4  2  5 )
          ( 5  4  3 )
     and the 5-by-2 matrix

          ( -10 -3 )
          (  12 14 )
      B = (  14 12 )
          (  16 16 )
          (  18 16 )
     We will first store the input matrix as a static C two-dimensional array,
     which is stored in row-major layout, and let LAPACKE handle the work space
     array allocation. The LAPACK base name for this function is gels, and we
     will use double precision (d), so the LAPACKE function name is LAPACKE_dgels.

     thus lda=3 and ldb=2. The output for each right hand side is stored in b as
     consecutive vectors of length 3. The correct answer for this problem is
     the 3-by-2 matrix

          ( 2 1 )
          ( 1 1 )
          ( 1 2 )
'''
print("\n(1b) - DGELS")
print("Solve the linear system ac = b, with\n")
print(">>> a = [[1,1,1],\n"
    + "         [2,3,4],\n"
    + "         [3,5,2],\n"
    + "         [4,2,5],\n"
    + "         [5,4,3]]")
print(">>> b = [[-10,3],\n"
    + "         [12,14],\n"
    + "         [14,12],\n"
    + "         [16,16],\n"
    + "         [18,16]]")
print(">>> A = Matrix('A',a)")
print(">>> B = Matrix('B',b)")
print(">>> C = et.dgels(A,B)")
print(">>> print(C)")
a = [[1,1,1],
     [2,3,4],
     [3,5,2],
     [4,2,5],
     [5,4,3]]
b = [[-10,3],
     [12,14],
     [14,12],
     [16,16],
     [18,16]]
A = Matrix('A',a)
B = Matrix('B',b)
C = et.dgels(A,B)
print(A)
print(B)
print(C)

'''
    Same example as (1b) except with Au = v where
    'u' and 'v' are vectors.
'''
print("\n(1c) - DGELS")
print("Solve the linear system Au = v, with\n")
print(">>> a = [[1,1,1],\n"
    + "         [2,3,4],\n"
    + "         [3,5,2],\n"
    + "         [4,2,5],\n"
    + "         [5,4,3]]")
print(">>> b = [1,2,3,4,5]")
print(">>> A = Matrix('A',a)")
print(">>> v = Matrix('v',b)")
print(">>> u = et.dgels(A,B)")
print(">>> print(C)")
a = [[1,1,1],
     [2,3,4],
     [3,5,2],
     [4,2,5],
     [5,4,3]]
b = [1,2,3,4,5]
A = Matrix('A',a)
v = Vector('v',b)
u = et.dgels(A,v)
print(A)
print(v)
print(u)

'''
    Same examples as (1b) and (1c) except using
    DGELSY as the driver.
'''
print("\n(1d) - DGELSY")
print("Solve the linear system ac = b, with\n")
print(">>> a = [[1,1,1],\n"
    + "         [2,3,4],\n"
    + "         [3,5,2],\n"
    + "         [4,2,5],\n"
    + "         [5,4,3]]")
print(">>> b = [[-10,3],\n"
    + "         [12,14],\n"
    + "         [14,12],\n"
    + "         [16,16],\n"
    + "         [18,16]]")
print(">>> A = Matrix('A',a)")
print(">>> B = Matrix('B',b)")
print(">>> C = et.dgelsy(A,B)")
print(">>> print(C)")
print(">>> rank_A = A.get_rank()")
print(">>> print('The rank of A is: ',rank_A)")
a = [[1,1,1],
     [2,3,4],
     [3,5,2],
     [4,2,5],
     [5,4,3]]
b = [[-10,3],
     [12,14],
     [14,12],
     [16,16],
     [18,16]]
A = Matrix('A',a)
B = Matrix('B',b)
C = et.dgelsy(A,B)
print(A)
print(B)
print(C)
rank_A = A.get_rank()
print('The rank of A is: ',rank_A)

print("\n(1e) - DGELSY")
print("Solve the linear system Au = v, with\n")
print(">>> a = [[1,1,1],\n"
    + "         [2,3,4],\n"
    + "         [3,5,2],\n"
    + "         [4,2,5],\n"
    + "         [5,4,3]]")
print(">>> b = [1,2,3,4,5]")
print(">>> A = Matrix('A',a)")
print(">>> v = Matrix('v',b)")
print(">>> u = et.dgelsy(A,B)")
print(">>> print(C)")
print(">>> rank_A = A.get_rank()")
print(">>> print('The rank of A is: ',rank_A)")
a = [[1,1,1],
     [2,3,4],
     [3,5,2],
     [4,2,5],
     [5,4,3]]
b = [1,2,3,4,5]
A = Matrix('A',a)
v = Vector('v',b)
u = et.dgelsy(A,v)
print(A)
print(v)
print(u)
rank_A = A.get_rank()
print('The rank of A is: ',rank_A)

'''
    Same examples as (1b) and (1c) except using
    DGELSD as the driver.
'''
print("\n(1f) - DGELSD")
print("Solve the linear system ac = b, with\n")
print(">>> a = [[1,1,1],\n"
    + "         [2,3,4],\n"
    + "         [3,5,2],\n"
    + "         [4,2,5],\n"
    + "         [5,4,3]]")
print(">>> b = [[-10,3],\n"
    + "         [12,14],\n"
    + "         [14,12],\n"
    + "         [16,16],\n"
    + "         [18,16]]")
print(">>> A = Matrix('A',a)")
print(">>> B = Matrix('B',b)")
print(">>> C = et.dgelsd(A,B)")
print(">>> print(C)")
print(">>> rank_A = A.get_rank()")
print(">>> sing_A = A.get_singular_values()")
print(">>> print('The rank of A is: ',rank_A)")
print(">>> print('The singular values of A are:\n',sing_A)")
a = [[1,1,1],
     [2,3,4],
     [3,5,2],
     [4,2,5],
     [5,4,3]]
b = [[-10,3],
     [12,14],
     [14,12],
     [16,16],
     [18,16]]
A = Matrix('A',a)
B = Matrix('B',b)
C = et.dgelsd(A,B)
print(A)
print(B)
print(C)
rank_A = A.get_rank()
sing_A = A.get_singular_values()
print('The rank of A is: ',rank_A)
print('The singular values of A are:\n',sing_A)

print("\n(1g) - DGELSD")
print("Solve the linear system Au = v, with\n")
print(">>> a = [[1,1,1],\n"
    + "         [2,3,4],\n"
    + "         [3,5,2],\n"
    + "         [4,2,5],\n"
    + "         [5,4,3]]")
print(">>> b = [1,2,3,4,5]")
print(">>> A = Matrix('A',a)")
print(">>> v = Matrix('v',b)")
print(">>> u = et.dgelsd(A,B)")
print(">>> print(C)")
print(">>> rank_A = A.get_rank()")
print(">>> print('The rank of A is: ',rank_A)")
a = [[1,1,1],
     [2,3,4],
     [3,5,2],
     [4,2,5],
     [5,4,3]]
b = [1,2,3,4,5]
A = Matrix('A',a)
v = Vector('v',b)
u = et.dgelsd(A,v)
print(A)
print(v)
print(u)
rank_A = A.get_rank()
sing_A = A.get_singular_values()
print('The rank of A is: ',rank_A)
print('The singular values of A are:\n',sing_A)

'''
    Same examples as (1b) and (1c) except using
    DGELSS as the driver.
'''
print("\n(1h) - DGELSS")
print("Solve the linear system ac = b, with\n")
print(">>> a = [[1,1,1],\n"
    + "         [2,3,4],\n"
    + "         [3,5,2],\n"
    + "         [4,2,5],\n"
    + "         [5,4,3]]")
print(">>> b = [[-10,3],\n"
    + "         [12,14],\n"
    + "         [14,12],\n"
    + "         [16,16],\n"
    + "         [18,16]]")
print(">>> A = Matrix('A',a)")
print(">>> B = Matrix('B',b)")
print(">>> C = et.dgelss(A,B)")
print(">>> print(C)")
print(">>> rank_A = A.get_rank()")
print(">>> sing_A = A.get_singular_values()")
print(">>> print('The rank of A is: ',rank_A)")
print(">>> print('The singular values of A are:\n',sing_A)")
a = [[1,1,1],
     [2,3,4],
     [3,5,2],
     [4,2,5],
     [5,4,3]]
b = [[-10,3],
     [12,14],
     [14,12],
     [16,16],
     [18,16]]
A = Matrix('A',a)
B = Matrix('B',b)
C = et.dgelss(A,B)
print(A)
print(B)
print(C)
rank_A = A.get_rank()
sing_A = A.get_singular_values()
print('The rank of A is: ',rank_A)
print('The singular values of A are:\n',sing_A)

print("\n(1i) - DGELSS")
print("Solve the linear system Au = v, with\n")
print(">>> a = [[1,1,1],\n"
    + "         [2,3,4],\n"
    + "         [3,5,2],\n"
    + "         [4,2,5],\n"
    + "         [5,4,3]]")
print(">>> b = [1,2,3,4,5]")
print(">>> A = Matrix('A',a)")
print(">>> v = Matrix('v',b)")
print(">>> u = et.dgelss(A,B)")
print(">>> print(C)")
print(">>> rank_A = A.get_rank()")
print(">>> print('The rank of A is: ',rank_A)")
a = [[1,1,1],
     [2,3,4],
     [3,5,2],
     [4,2,5],
     [5,4,3]]
b = [1,2,3,4,5]
A = Matrix('A',a)
v = Vector('v',b)
u = et.dgelss(A,v)
print(A)
print(v)
print(u)
rank_A = A.get_rank()
sing_A = A.get_singular_values()
print('The rank of A is: ',rank_A)
print('The singular values of A are:\n',sing_A)

print("\n(2) DGETRF")

print("\n(2a) - DGETRF")
print(">>> A = Matrix('A',[[1,2,3],[4,5,6],[7,8,9]])")
print(">>> pivots = et.dgetrf(A)")
print(">>> perm = et.permutation_matrix(A.get_num_rows,pivots)")
print(">>> print(pivots)")
print(">>> print(perm)")
A = Matrix('A',[[1,2,3],[4,5,6],[7,8,9]])
pivots = et.dgetrf(A)
perm = et.permutation_matrix(A.get_num_rows(),pivots)
print(pivots)
print(perm)

print("\n(2b) - DGETRF_LU")
print(">>> A = Matrix('A',[[1,2,3],[4,5,6],[7,8,9]])")
print(">>> LU = et.dgetrf_l_u(A)")
print(">>> print(LU)")
A = Matrix('A',[[1,2,3],[4,5,6],[7,8,9]])
LU = et.dgetrf_l_u(A)
print(LU)

print("\n(2c) - DGETRF_L_U")
print(">>> A = Matrix('A',[[1,2,3],[4,5,6],[7,8,9]])")
print(">>> L,U = et.dgetrf_lu(A)")
print(">>> LU = L*U")
print(">>> print(A)")
print(">>> print(L)")
print(">>> print(U)")
print(">>> print(LU)")
A = Matrix('A',[[1,2,3],[4,5,6],[7,8,9]])
L,U = et.dgetrf_lu(A)
LU = L*U
print(A)
print(L)
print(U)
print(LU)

print("\n(2d) - DGETRF_PLU")
print(">>> A = Matrix('A',[[1,2,3],[4,5,6],[7,8,9]])")
print(">>> P,L,U = et.dgetrf_plu(A)")
print(">>> PLU = P*L*U")
print(">>> print(A)")
print(">>> print(L)")
print(">>> print(U)")
print(">>> print(LU)")
A = Matrix('A',[[1,2,3],[4,5,6],[7,8,9]])
P,L,U = et.dgetrf_plu(A)
PLU = P*L*U
print(A)
print(L)
print(U)
print(PLU)

print("\n(3) DGEQRF")
print("\n(3a) - DGEQRF")
print(">>> A = Matrix('A',[[1,2,3],[4,5,6],[7,8,9]])")
print(">>> ref = et.dgeqrf(A)")
print(">>> print(ref)")
A = Matrix('A',[[1,2,3],[4,5,6],[7,8,9]])
ref = et.dgeqrf(A)
print(ref)

print("\n(3b) - DORGQR")
print(">>> A = Matrix('A',[[1,2,3],[4,5,6],[7,8,9]])")
print(">>> ref = et.dgeqrf(A)")
print(">>> Q = et.dorgqr(A,ref)")
print(">>> print(ref)")
print(">>> print(Q)")
A = Matrix('A',[[1,2,3],[4,5,6],[7,8,9]])
ref = et.dgeqrf(A)
Q = et.dorgqr(A,ref)
print(ref)
print(Q)

print("\n(3c) - DGEQRF_QR")
print(">>> A = Matrix('A',[[1,2,3],[4,5,6],[7,8,9]])")
print(">>> Q,R = et.dgeqrf_qr(A)")
print(">>> QR = Q*R")
print(">>> print(A)")
print(">>> print(Q)")
print(">>> print(R)")
print(">>> print(QR)")
A = Matrix('A',[[1,2,3],[4,5,6],[7,8,9]])
Q,R = et.dgeqrf_qr(A)
QR = Q*R
print(A)
print(Q)
print(R)
print(QR)

'''
    Example from https://www.ibm.com/support/knowledgecenter/SSFHY8_6.2/reference/am5gr_hgesvd.html#am5gr_hgesvd__gesvddt

            |  1.0    1.0    0.0    0.0 |
    A    =  |  0.0    2.0    1.0    0.0 |
            |  0.0    0.0    3.0    1.0 |
            |  0.0    0.0    0.0    4.0 |

    Copy code

    Output:

    Array A is overwritten.


            | 4.260007 |
    S    =  | 3.107349 |
            | 2.111785 |
            | 0.858542 |


    INFO = 0
'''
print("\n(4) DGESVD")
print("\n(4a) - DGESVD")
print(">>> a = [[1, 1, 0, 0],\n"
    + "         [0, 2, 1, 0],\n"
    + "         [0, 0, 3, 1],\n"
    + "         [0, 0, 0, 4]]")
print(">>> A = Matrix('A',a)")
print(">>> s = et.dgesvd(A)")
print(">>> print(s)")
a = [[1,1,0,0],[0,2,1,0],[0,0,3,1],[0,0,0,4]]
A = Matrix('A',a)
s = et.dgesvd(A)
print(s)

print("\n(4b) - DGESDD")
print(">>> a = [[1, 1, 0, 0],\n"
    + "         [0, 2, 1, 0],\n"
    + "         [0, 0, 3, 1],\n"
    + "         [0, 0, 0, 4]]")
print(">>> A = Matrix('A',a)")
print(">>> s = et.dgesdd(A)")
print(">>> print(s)")
a = [[1,1,0,0],[0,2,1,0],[0,0,3,1],[0,0,0,4]]
A = Matrix('A',a)
s = et.dgesvd(A)
print(s)

'''
    Example fromhttps://www.ibm.com/support/knowledgecenter/SSFHY8_6.2/reference/am5gr_hgesvd.html#am5gr_hgesvd__gesvddt

            |  1.0    2.0    3.0  |
    A    =  |  2.0    4.0    5.0  |
            |  3.0    5.0    6.0  |

    Copy code

    Output:

    Array A is overwritten.


            | 11.344814 |
    S    =  |  0.515729 |
            |  0.170915 |

    Copy code


            | -0.327985 -0.736976 -0.591009 |
    U    =  | -0.591009 -0.327985  0.736976 |
            | -0.736976  0.591009 -0.327985 |

    Copy code


            | -0.327985 -0.591009 -0.736976 |
    VT   =  |  0.736976  0.327985 -0.591009 |
            | -0.591009  0.736976 -0.327985 |


    INFO = 0
'''
print("\n(4c) - DGESVD_SVD")
print(">>> a = [[1, 2, 3],\n"
    + "         [2, 4, 5],\n"
    + "         [3, 5, 6]]")
print(">>> A = Matrix('A',a)")
print(">>> U, S, VT = et.dgesvd_svd(A)")
print(">>> print(A)")
print(">>> print(U)")
print(">>> print(S)")
print(">>> print(VT)")
print(">>> B = U * S * VT")
print(">>> print(B)")
a = [[1, 2, 3],[2, 4, 5],[3, 5, 6]]
A = Matrix('A',a)
U, S, VT = et.dgesvd_svd(A)
print(A)
print(U)
print(S)
print(VT)
B = U * S * VT
print(B)

print("\n(4d) - DGESDD_SVD")
print(">>> a = [[1, 2, 3],\n"
    + "         [2, 4, 5],\n"
    + "         [3, 5, 6]]")
print(">>> A = Matrix('A',a)")
print(">>> U, S, VT = et.dgesdd_svd(A)")
print(">>> print(A)")
print(">>> print(U)")
print(">>> print(S)")
print(">>> print(VT)")
print(">>> B = U * S * VT")
print(">>> print(B)")
a = [[1, 2, 3],[2, 4, 5],[3, 5, 6]]
A = Matrix('A',a)
U, S, VT = et.dgesvd_svd(A)
print(A)
print(U)
print(S)
print(VT)
B = U * S * VT
print(B)

'''
    Example from

     A    = |  1.0    2.0    3.0    4.0  |
            |  5.0    6.0    7.0    8.0  |


    Output:

        For DGESVD:


        A    =  | -0.352062 -0.443626 -0.535190 -0.626754 |
                |  0.758981  0.321242 -0.116498 -0.554238 |

        For DGESDD, A has been overwritten on output.


    S    =  | 14.227407 |
            |  1.257330 |


    U    =  | -0.376168 -0.926551 |
            | -0.926551  0.376168 |

        For DGESVD, VT is not referenced.
        For DGESDD, VT is:


        VT   =  | -0.352062 -0.443626 -0.535190 -0.626754 |
                |  0.758981  0.321242 -0.116498 -0.554238 |


    INFO = 0
'''
print("\n(4e) - DGESVD_SVD")
print(">>> a = [[1, 2, 3, 4],\n"
    + "         [5, 6, 7, 8]]")
print(">>> A = Matrix('A',a)")
print(">>> U, S, VT = et.dgesvd_svd(A)")
print(">>> print(A)")
print(">>> print(U)")
print(">>> print(S)")
print(">>> print(VT)")
print(">>> B = U * S * VT")
print(">>> print(B)")
a = [[1, 2, 3, 4],
     [5, 6, 7, 8]]
A = Matrix('A',a)
U, S, VT = et.dgesvd_svd(A)
print(A)
print(U)
print(S)
print(VT)
B = U * S * VT
print(B)

print("\n(4f) - DGESDD_SVD")
print(">>> a = [[1, 2, 3, 4],\n"
    + "         [5, 6, 7, 8]]")
print(">>> A = Matrix('A',a)")
print(">>> U, S, VT = et.dgesdd_svd(A)")
print(">>> print(A)")
print(">>> print(U)")
print(">>> print(S)")
print(">>> print(VT)")
print(">>> B = U * S * VT")
print(">>> print(B)")
a = [[1, 2, 3, 4],
     [5, 6, 7, 8]]
A = Matrix('A',a)
U, S, VT = et.dgesvd_svd(A)
print(A)
print(U)
print(S)
print(VT)
B = U * S * VT
print(B)
