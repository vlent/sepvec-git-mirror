function C = sepvec_diff(A, B)

C.coeff = [ A.coeff, -B.coeff ];
C.vec = [ A.vec, B.vec ];
