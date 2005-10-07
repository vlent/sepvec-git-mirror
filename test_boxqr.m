% number of dimensions
d = 3;
% size of vectors (same for each dimension)
M = 100;
% seperation ranks
r = 6;
s = r + 3;
for i = 1:d
  A{i} = rand(M, r);
  B{i} = rand(M, s);
end

% (1) the products are multiplied out 
% -> exponential in d
tic
AA = A{1};
BB = B{1};
for i = 2:d
  AA = boxprod(AA, A{i});  
  BB = boxprod(BB, B{i});
end

[CC, RR] = qraug(AA, BB);
toc

% (2) the products are eliminated pairwise
tic
[C_, R_] = boxqr_(A, B);
toc

% (3) QR on matrices first, then eliminate pairwise
tic
[C, R] = boxqr(A, B);
toc

% difference in the least-squares solutions
err = [maxabs(RR\CC - R_\C_)
       maxabs(RR\CC - R\C)
       maxabs(R_\C_ - R\C)]';
digits = -log10(err)