function Y = invveck(X,p)
k = 0;
Y = zeros(p,p);

for j = 1 : p-1
    for i = j+1:p
         k = k + 1;
         Y(i,j)=X(k);
    end
end
    
Y = Y-Y';

end