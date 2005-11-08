function [dims, r, s] = sepvec_compat(F, G)

dims = length(F.coeff);
if dims != length(G.coeff), error('number of dimensions mest be equal'); end
r = size(F.vec{1}, 2);
s = size(G.vec{1}, 2);
for i = 2:dims
  if r != size(F.vec{i}, 2), error('seperation rank of first argument'); end
  if s != size(G.vec{i}, 2), error('seperation rank of second argument'); end
end
for i = 1:dims
  if size(F.vec{i}, 1) != size(G.vec{i}, 1)
    error('length of corresponding dimensions must be equal'); 
  end
end
