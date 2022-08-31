function [H11,H12,H22] = makeH(MZ,MZo,MZoo,N,n,p)
  Ip = eye(p);
  Inmp = eye(n-p);
  S11 = zeros(p^2,p^2);
  S12 = zeros(p^2,p*(n-p));
  H22 = zeros(p*(n-p),p*(n-p));
  Deltap = makeDelta(p);
  Dp = makeD(p);
  for l = 1 : N
      Zl = mget(MZ,l,p);
      Zlo = mget(MZo,l,n-p);
      Zloo = mget(MZoo,l,n-p);
      D = diag(diag(Zl));
      kIpZl = kron(Ip,Zl);
      S11 = S11 + kron(D,Zl) + 2*kIpZl*Deltap*kIpZl-kron(symm(Zl*D),Ip);
      S12 = S12 + kron(D,Zlo) + 2*kIpZl*Deltap*kron(Ip,Zlo);
      H22 = H22 + kron(D,Zloo)+2*kron(Ip,Zlo')*Deltap*kron(Ip,Zlo)-kron(symm(Zl*D),Inmp);
  end
 H11 = symm(-2*Dp'*S11*Dp);
 H12 = -2*Dp'*S12;
 H22 = symm(-4*H22);
end