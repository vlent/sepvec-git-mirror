function C = sepvec_sub(A, B)

C.coeff = [ A.coeff, -B.coeff ];
C.vec = [ A.vec, B.vec ];
