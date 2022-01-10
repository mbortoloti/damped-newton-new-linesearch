%
%   This function reads initial guesses from file given in 'address'
%   with a formmat given by 'formmat'
%
function [A] = getinput(address,formatSpec,n,m)
    fileID = fopen(address,'r');
    cellA = textscan(fileID,formatSpec,n*m);
    k = 0;
    for i = 1 : n
        for j = 1 : m
            k = k + 1;
            A(i,j) = cellA{1}(k);
        end
    end
    fclose(fileID);
end