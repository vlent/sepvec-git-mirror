function [C, R] = boxqr_(A, B)

d = length(A);
r = size(A{1}, 2);

% pairwise elimination of the `box products'
R = A{1};
C = B{1};
for i = 2:d
  AA = boxprod(R, A{i});
  BB = boxprod(C, B{i});
  [C, R] = qraug(AA, BB);
end
