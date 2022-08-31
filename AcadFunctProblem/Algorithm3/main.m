%
%   An efficient damped Newton-type algorithm with
%   globalization strategy on Riemannian manifolds
%
%   by Bortoloti, M. A. A., Fernandes, T. A. and Ferreira, O. P.
%
%   Algorithm 3 developed to minimize
%
%  f(P) = a log(det(P)) + b trace(P^-1),
%
%  where a,b > 0, on the  cone of symmetric positive definite matrices.
%   
%
clear all;
clc;


% Dimension setting
n = 5;
fprintf("n = %5d\n",n);

% linesearch parameter setting
theta = 0.1;

% Definition of Retraction 
%
%      (Choose one of the following options 
%       in order to define the retraction)
%
% typeretra == 1 Exponential map
% typeretra == 2 Analogous Exponential map
% typeretra == 3 Second order approximarion for exponential map
% typeretra == 4 First order approximarion for exponential map
typeretra = 1;


% Initial guess for testing algorithm 1
density = rand(1);
rc = rand(1);
P = sprandsym(n,density,rc,1);
P= full(P);


options.maxiter = 2000;
options.stpmin = 1.e-10;
options.ngtol = 1.0e-6;
options.eps2 = 1.0e-16;
 
 
options.a = 3.0;
options.b = 1.0;

        
I = eye(n);
        
% Riemannian gradient of f at P
a = options.a;
b = options.b;
rgradf = @(P) a*P-b*I;
options.rgradf = rgradf;
      
% Retraction definition
switch typeretra
        case 1 
             options.ret = @exponential;
        case 2
            options.ret = @exponential2;
        case 3
            options.ret = @secondorder;
        case 4
            options.ret = @firstorder;
        otherwise
            fprintf('Retraction not defined.\n');
end

        
% Metric definition (Rothaus metric)
metric = @(U,V,P) trace(V*P^(-1)*U*P^(-1));
options.metric = metric;
 
% Solver call 
[P,info] = algorithm3(P,options,theta);

%
% Retractions
%
function Y = firstorder(P,V)
    Y = symm(P+V);
end

function Y = exponential(P,V)
    Y = sqrtm(P)*expm(sqrtm(P)^(-1)*V*sqrtm(P)^(-1))*sqrtm(P); 
end

function Y = exponential2(P,V)
 Y = symm(P*expm(P^(-1)*V));
end

function Y = secondorder(P,V)
    Y = symm(P + V + 0.5*V*P^(-1)*V);
end
function Y = symm(X)
    Y = 0.5*(X+X');
end