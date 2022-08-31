function [H] = LHSnewton(MZ,MZo,MZoo,N,n,p)
 [H11,H12,H22] = makeH(MZ,MZo,MZoo,N,n,p);
 H = [H11,H12;2*H12',H22];
end