function [C, R] = boxqr_reg(A, B, mu)

% A, B are lists of matrices with the same numbers of columns
% the corresponding matrices in A and B have the same number of rows
% R is the lower triangular matrix in the QR-factorization of
% the box product of the matrices in A
% C is Q' applied to the `box product' of the matrices in B
%
% it is assumed that the matrices have no less rows than columns

d = length(A);
r = size(A{1}, 2);

% pairwise elimination of the `box products'
[C, R] = qraug_reg(A{1}, B{1}, mu);
for i = 2:d
  [C_, R_] = qraug_reg(A{i}, B{i}, mu);
  RR = boxprod(R, R_);
  CC = boxprod(C, C_);
  %[C, R] = qraug(RR, CC);
  [C, R] = qraug_reg(RR, CC, mu);
end
