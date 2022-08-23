function [V] = vec(X)
    [~,p] = size(X);
    V = X(:,1);
    for i = 2 : p
        V = [V;X(:,i)];
    end
end