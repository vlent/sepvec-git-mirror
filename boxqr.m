function [C, R] = boxqr(A, B)

% boxqr   QR-factorization of box product
%
% [C, R] = boxqr(A, B)
% A and B are lists of matrices so that for all i and j
%   size(A{i}, 2) == size(A{j}, 2)
%   size(B{i}, 2) == size(B{j}, 2)
%   size(A{i}, 1) == size(B{i}, 1)
% R is the lower triangular matrix in the reduced QR-factorization of
% the `box product' of the matrices in A
% C is Q' applied to the `box product' of the matrices in B
%
% it is assumed that the matrices have no less rows than columns
%   size(A{i}, 1) >= size(A{i}, 2)

d = length(A);
r = size(A{1}, 2);

% pairwise elimination of the `box products'
[C, R] = qraug(A{1}, B{1});
for i = 2:d
  [C_, R_] = qraug(A{i}, B{i});
  RR = boxprod(R, R_);
  CC = boxprod(C, C_);
  [C, R] = qraug(RR, CC);
end
