function Y = invvec(X,m1,m2)
  Y = X(1:m1);
  for j = 2 : m2
      k1 = (j-1)*m1+1;
      k2 = j*m1;
      Y = [Y,X(k1:k2)];
  end
end