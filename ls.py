from mynum import *

def norm(x):
    return maximum.reduce(absolute(x).flat)

m = 100
n = 40
A_ = uniform(-1.0, 1.0, (m, n))
U, s_, Vt = singular_value_decomposition(A_)
print norm(matrixmultiply(matrixmultiply(U, identity(n) * s_), Vt) - A_)
s[0] *= 1000
A =matrixmultiply(matrixmultiply(U, identity(n) * s), Vt)
#x = uniform(-1.0, 1.0, n)
b = uniform(-1.0, 1.0, m)

x1, r, rank, s1 = linear_least_squares(A, b)
At = transpose(A)
x2 = solve_linear_equations(matrixmultiply(At, A), matrixmultiply(At, b))
print norm(x1-x2)
