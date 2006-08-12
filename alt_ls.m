function [F, it, flag, norm_res] = alt_ls(G, F, maxit, tol, mu)

%%% F = alt_ls(G, F)
%%% 
%%% alternating least squares algorithm
%%% given a vector G in separated form 
%%% (list of matrices with same number of columns)
%%% find an approximation F in separated form
%%% starting from given F
%%% maxit: maximum number of iterations
%%% tol: relative tolerance
%%% mu: regularization parameter
%%% it: number of iterations performed
%%% flag: 1 when succesful, 0 when maxit reached before convergence
%%% norm_res: vector with norms of residuals

[dims, r, s] = sepvec_compat(F, G);

nG = sepvec_norm(G);
it = 0;
flag = 0;
norm_res(it+1) = sepvec_norm(sepvec_sub(G, F));
%%% this whole loop can be more efficient by saving intermediate results
while (!flag)
  if norm_res(it+1) <= tol * nG
    flag = 1
  elseif (it < maxit)
    for i = 1:dims
      F = alt_ls_step(G, F, i, mu);
    end
    it = it + 1;
    norm_res(it+1) = sepvec_norm(sepvec_sub(G, F));
  else
    break;
  end
end

