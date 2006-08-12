dims = 3;
r = 4;
s = 4;
M = [ 5 7 5 ];

F.coeff = rand(r, 1);
G.coeff = rand(s, 1);
F.vec = {};
G.vec = {};
for i = 1:dims
  F.vec{i} = rand(M(i), r);
  F.vec{i} = F.vec{i} / diag(sqrt(sum(F.vec{i}.^2, 1)));
  G.vec{i} = rand(M(i), s);
  G.vec{i} = G.vec{i} / diag(sqrt(sum(G.vec{i}.^2, 1)));
end
[ dims_, r_, s_ ] = sepvec_compat(F, G);
all([dims, r, s] == [dims_, r_, s_])
H = sepvec_sub(F, G);
sepvec_norm(F)
sepvec_norm(sepvec_sub(F, F))
sepvec_norm(sepvec_sub(G, G))

%F2.coeff = rand(r, 1);
pert = 1e-2;
F2.coeff = F.coeff + pert * randn(r, 1);
F2.vec = {};
for i = 1:dims
  %F2.vec{i} = rand(M(i), r);
  F2.vec{i} = F.vec{i} + pert * randn(M(i), r);
  F2.vec{i} = F2.vec{i} / diag(sqrt(sum(F2.vec{i}.^2, 1)));
end

maxit = 10;
tol = 1e-6;
mu = 0;
sepvec_norm(sepvec_sub(F, F2))
[F3, it, flag, norm_res] = alt_ls(F, F2, maxit, tol, mu);
it
flag
sepvec_norm(sepvec_sub(F, F3))
