function [C, R] = qraug_reg(A, B, mu)

% qraug_reg   regularized reduced QR-factorization of augmented matrix
% 
% [C, R] = qraug_reg(A, B, mu)

r = size(A, 2);
s = size(B, 2);
[C, R] = qraug([A; eye(r)*mu], [B; zeros(r, s)])



