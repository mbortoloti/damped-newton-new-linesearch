function X = comporth(X)
    [m,p]=size(X);
    [XX,~] = qr(X);
    X = XX(:,p+1:m);
end