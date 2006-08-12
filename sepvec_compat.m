function [dims, r, s] = sepvec_compat(F, G)

dims = length(F.vec);
if dims != length(G.vec), error('number of dimensions must be equal'); end
r = size(F.coeff, 1);
s = size(G.coeff, 1);
for i = 1:dims
  if r != size(F.vec{i}, 2), error('separation rank of first argument'); end
  if s != size(G.vec{i}, 2), error('separation rank of second argument'); end
end
for i = 1:dims
  if size(F.vec{i}, 1) != size(G.vec{i}, 1)
    error('length of corresponding dimensions must be equal'); 
  end
end
