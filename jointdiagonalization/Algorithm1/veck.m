function [V] = veck(X)
[~,p] = size(X);
    k = 0;
    for j = 1 : p
        for i = j+1 : p
            k = k + 1;
            V(k,1) = X(i,j);
        end
    end
end
    