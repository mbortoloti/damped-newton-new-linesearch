function [X] = mget(MX,l,ncol)
  X = MX(:,(l-1)*ncol+1:l*ncol);
end