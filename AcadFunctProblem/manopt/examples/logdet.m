%
%   logdet problem
%


%
%  objective function
%
a = 1.0;
b = 1.0;
f = @(X) a * log(det(X)) + b * trace(X^-1);


