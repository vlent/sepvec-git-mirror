import math

from scipy import *
from linalg import *


def boxprod(A, B):
    mA, nA = A.shape
    mB, nB = B.shape
    assert nA == nB
    return reshape(A[:,NewAxis,:] * B[NewAxis,:,:], (mA * mB, nA))

class symprod:
    def __init__(self, name):
        self.name = name
    def __mul__(self, other):
        return symprod(self.name + "*" + other.name)
    def __str__(self):
        return self.name
    def __repr__(self):
        return self.name

def test_boxprod():
    A = reshape(map(symprod, map(str, arange(3*2)+1)),(3,2))
    B = reshape(map(symprod, map(str, arange(2*2)+1)),(2,2))
    C = boxprod(A, B)
    print C

class sepvec:
    def __init__(self, coeff, vec):
        assert len(coeff) == len(vec)
        self.seprank = len(coeff)
        assert alltrue([ v.shape[-1] == self.seprank for v in vec ])
        self.dim = [ v.shape[:-1] for v in vec ]
        self.coeff = coeff
        self.vec = vec

    def __add__(self, other):
        return sepvec(self.coeff + other.coeff, self.vec + other.vec)

    def __neg__(self):
        neg_coeff = [ -c for c in self.coeff ]
        return sepvec(neg_coeff, self.vec)

    def __sub__(self, other):
        return self + (-other)

    def __dot__(self, other):
        n, r, s = self.compat(other)
        X = 1
        for i in range(n):
            X = X * dot(self.vec[i], self.vec[i])
        return dot(dot(self.coeff, X), other.coeff)

    def __norm__(self):
        return math.sqrt(abs(self.dot(self)))

def qr_r(A):
    lwork = qr(A, overwrite_a=1, lwork=-1)
    # builds Q, not necessary, try to find a reduced QR algorithm
    Q, R = qr(A, overwrite_a=1, lwork)
    return R

def qraug_concat(A, B):
    AB = concatenate((A, B), -1)
    RC = qr_r(AB)
    return RC

def qraug(A, B):
    RC = qraug_concat(A, B)
    r = A.shape[-1]
    return RC[:,:r], RC[:,r:]

def boxqr(A, B):
    R, C = qraug(A[0], B[0])
    for i in range(1, len(A)):
        R_, C_ = qraug(A[i], B[i])
        RR = boxprod(R, R_)
        CC = boxprod(C, C_)
        R, C = qraug(RR, CC)
    return R, C

def boxqr_(A, B):
    RC0 = qraug_concat(A[0], B[0])
    for i in range(1, len(A)):
        RC1 = qraug_concat(A[i], B[i])
        RC2 = boxprod(RC0, RC1)
        RC0 = qr_r(RC2)
    r = A.shape[-1]
    R = RC0[:,:r]
    C = RC0[:,r:]
    return R, C
