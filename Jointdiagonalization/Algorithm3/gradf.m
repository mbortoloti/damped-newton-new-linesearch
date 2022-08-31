function [Gf] = gradf(Y,MA,N,n,p)
 Gf = zeros(n,p);
 for l = 1 : N
    Al = mget(MA,l,n);
    D = diag(diag(Y'*Al*Y));
    Gf = Gf + Al*Y*D;
 end
 Gf = -4*Gf;
 Gf = proj(Y,Gf);
end