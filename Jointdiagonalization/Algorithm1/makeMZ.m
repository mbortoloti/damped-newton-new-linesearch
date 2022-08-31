function [MZ] = makeMZ(MA,Y1,Y2,N,n)
%  Zl = Y1'*Al*Y2
Al = mget(MA,1,n);
MZ = Y1'*Al*Y2;
for l = 2:N
  Al = mget(MA,l,n);  
  MZ = [MZ,Y1'*Al*Y2];
end


end