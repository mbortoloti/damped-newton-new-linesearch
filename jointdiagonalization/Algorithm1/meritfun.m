function [r] = meritfun(Y,MA,N,n,p)
  G = gradf(Y,MA,N,n,p);
  r = 0.5*norm(G,'fro')^2;
end