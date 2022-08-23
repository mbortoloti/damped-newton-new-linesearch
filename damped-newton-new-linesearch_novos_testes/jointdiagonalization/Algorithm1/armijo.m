function [alpha,evalf] = armijo(Y,V,Gm,MA,N,n,p,sigma,evalf)
  alpha = 1;
  phiY = meritfun(Y,MA,N,n,p);
  evalf = evalf + 1;
  inndot = trace(Gm'*V);
  while 1
      Q = qf(Y+alpha*V);
      phiQ = meritfun(Q,MA,N,n,p);
      evalf = evalf + 1;
      arm = phiQ-phiY-sigma*alpha*inndot;
      if arm > 0
          alpha = 0.5*alpha;
      else
          break;
      end
  end
end