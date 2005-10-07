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

% QR-factorization of the individual matrices
%for i = 1:d
%  [CC{i}, RR{i}] = qr(A{i}, B{i}, 0);
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
