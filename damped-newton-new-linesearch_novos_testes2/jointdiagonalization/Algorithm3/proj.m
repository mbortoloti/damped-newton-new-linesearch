function Z = proj(X,Y)
  % Z = P_X(Y) = Y-Xsym(X'*Y)
  Z = Y - X*symm(X'*Y);
end
