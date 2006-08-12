function C = sepvec_sub(A, B)

[dims, r, s] = sepvec_compat(A, B);
C.coeff = [ A.coeff; -B.coeff ];
for i = 1:dims
  C.vec{i} = [ A.vec{i}, B.vec{i} ];
end
