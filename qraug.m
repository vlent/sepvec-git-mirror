function [C, R] = qraug(A, B)

% qraug   reduced QR-factorization of augmented matrix
% 
% [C, R] = qraug(A, B)

RC = qr([A, B], 0);
n = size(A, 2);
R = triu(RC(1:n, 1:n));
C = RC(1:n, n+1:end);


