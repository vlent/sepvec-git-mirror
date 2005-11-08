function [F, it, flag] = alt_ls(G, F, maxit, tol, mu)

%%% F = alt_ls(G, F)
%%% 
%%% alternating least squares algorithm
%%% given a vector G in seperated form 
%%% (list of matrices with same number of columns)
%%% find an approximation F in seperated form
%%% starting from given F
%%% maxit: maximum number of iterations
%%% tol: relative tolerance
%%% mu: regularization parameter
%%% it: number of iterations performed
%%% 

[dims, r, s] = sepvec_compat(F, G);

nG = sepvec_norm(G);
it = 0;
flag = 0;
%%% this whole loop can be more efficient by saving intermediate results
while (!flag) && (it < maxit)
  if sepvec_norm(sep_vec_sub(G, F)) <= tol * nG
    flag = 1
  else
    for i = 1:dims
      F = alt_ls_step(G, F, i, mu);
    end
    it = it + 1;
  end
end

