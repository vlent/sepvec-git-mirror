F = alt_ls_step(G, F, i, mu)

[dims, r, s] = sepvec_compat(F, G);

C = ones(1, r);
R = ones(1, r);
for j = 1:dims
  if j != i
    [C_, R_] = qraug(F.vec{j}, G.vec{j});
    RR = boxprod(R, R_);
    CC = boxprod(C, C_);
    [C, R] = qraug(RR, CC);
  end
end
if mu != 0
  [C, R] = qraug([R; mu*eye(r)], [C; zeros(r, s)]);
end

% boxprod(R, F.vec{i}) * F.coeff == boxprod(C, G.vec{i}) * G.coeff
% sum(boxprod(R, F.vec{i} * diag(F.coeff)), 2)
%   == sum(boxprod(C, G.vec{i} * diag(G.coeff)), 2)
% R * (F.vec{i} * diag(F.coeff))' == C * (G.vec{i} * diag(G.coeff))'
% R * X' == C * (G.vec{i} * diag(G.coeff))'

Xt = (R \ C) * (G.vec{i} * diag(G.coeff))';
%X = G.vec{i} * ((diag(G.coeff) * C') / R');
c = sqrt(sum(Xt.^2, 2));
%F.vec{i} = X / diag(c); 
F.vec{i} = dmult(1 ./ c, Xt)';
F.coeff = c;
