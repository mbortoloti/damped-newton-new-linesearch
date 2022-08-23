function [D] = makeDelta(n)
D = zeros(n^2,n^2);
for i = 1 : n
   D = D + kron(En(i,i,n),En(i,i,n)); 
end


    function [E] = En(i,j,n)
        E = zeros(n,n);
        E(i,j) = 1.0;
    end

end