function C = boxprod(A, B)

% return the `box product' of two matrices
% the resulting matrix has as columns the tensor products of 
% the corresponding columns of the original matrices
% A \boxtimes B = [ a_1 \otimes b_1, \dotsc, a_r \otimes b_r ]


C = kron(A, ones(size(B, 1), 1)) .* kron(ones(size(A, 1), 1), B);
