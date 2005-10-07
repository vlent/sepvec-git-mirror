function [C, R] = boxqr_(A, B)

% A, B are lists of matrices with the same numbers of columns
% the corresponding matrices in A and B have the same number of rows
% R is the lower triangular matrix in the QR-factorization of
% the box product of the matrices in A
% C is Q' applied to the `box product' of the matrices in B

d = length(A);
r = size(A{1}, 2);

% pairwise elimination of the `box products'
R = A{1};
C = B{1};
for i = 2:d
  AA = boxprod(R, A{i});
  BB = boxprod(C, B{i});
  RC = qr([AA, BB], 0);
  R = triu(RC(1:r,1:r));
  C = RC(1:r,r+1:2*r);
end
