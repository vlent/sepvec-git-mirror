function [C, R] = boxqr(A, B)

% A, B are lists of matrices with the same numbers of columns
% the corresponding matrices in A and B have the same number of rows
% R is the lower triangular matrix in the QR-factorization of
% the box product of the matrices in A
% C is Q' applied to the `box product' of the matrices in B
%
% it is assumed that the matrices have no less rows than columns

d = length(A);
r = size(A{1}, 2);

% QR-factorization of the individual matrices
%for i = 1:d
%  [CC{i}, RR{i}] = qr(A{i}, B{i});
%end

% pairwise elimination of the `box products'
RC = qr([A{1}, B{1}], 0);
R = triu(RC(1:r,1:r));
C = RC(1:r,r+1:2*r);
for i = 2:d
  RC_ = qr([A{i}, B{i}], 0);
  R_ = triu(RC_(1:r,1:r));
  C_ = RC_(1:r,r+1:2*r);
  RR = boxprod(R, R_);
  CC = boxprod(C, C_);
  RC = qr([RR, CC], 0);
  R = triu(RC(1:r,1:r));
  C = RC(1:r,r+1:2*r);
end
