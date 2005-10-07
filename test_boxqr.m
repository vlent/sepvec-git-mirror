d = 3;
M = 100;
r = 7;
for i = 1:d
  A{i} = rand(M, r);
  B{i} = rand(M, r);
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

RRCC = qr([AA, BB], 0);
RR = triu(RRCC(1:r,1:r));
CC = RRCC(1:r,r+1:2*r);
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
RR\CC - R_\C_
RR\CC - R\C
