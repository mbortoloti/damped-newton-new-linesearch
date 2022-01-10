function [Gm] = gradmerit(Y,Gf,MA,N,n,p)
  Gm = zeros(n,p);
  for l = 1 : N
      Al = mget(MA,l,n);
      DYtAlY = diag(diag(Y'*Al*Y));
      Gm = Gm + Al*Gf*DYtAlY+2*Al*Y*diag(diag(Y'*Al*Gf)) ...
              -Gf*symm(Y'*Al*Y*DYtAlY);
  end
  Gm = -4*proj(Y,Gm);
end