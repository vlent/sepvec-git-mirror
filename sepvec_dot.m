function a = sepvec_dot(F, G)

[dims, r, s] = sepvec_compat(F, G);

%X = ones(r, s);
X = 1;
for i = 1:dims
  X = X .* (F.vec{i}' * G.vec{i});
end
a = F.coeff' * X * G.coeff;
