function [G] = RHSnewton(MZ,MZo,N,n,p)
G1 = zeros(p,p);
G2 = zeros(n-p,p);
  for l = 1:N
      Zl = mget(MZ,l,p);
      D = diag(diag(Zl));
      G1 = G1 + Zl*D;
      Zlo = mget(MZo,l,n-p);
      G2 = G2 + Zlo'*D;
  end
  G1 = -4*skew(G1);
  G2 = -4*G2;
  G = [veck(G1);vec(G2)];
end