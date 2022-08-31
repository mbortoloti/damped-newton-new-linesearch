function [D] = makeD(n)
    D = zeros(n^2,n*(n-1)/2);
    for i = 1 : n
        for j = 1 : i-1
          D = D + En(n*(j-1)+i,j*(n-(j+1)/2)-n+i,n^2,n*(n-1)/2) ...
                - En(n*(i-1)+j,j*(n-(j+1)/2)-n+i,n^2,n*(n-1)/2);            
        end
        
    end



    function [E] = En(i,j,p,q)
        E = zeros(p,q);
        E(i,j) = 1;
    end

end