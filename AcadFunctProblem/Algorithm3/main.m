%
%   On the globalization of Riemannian Newton method
%   by Bortoloti, M. A. A., Fernandes, T. A. and Ferreira, O. P.
%
%   Algorithm 3 developed to solve an academic problem on the 
%   cone of symmetric positive definite matrices truncated 
%   singular value for f_1 function 
%   
%
clear all;
clc;

% Dimension definition 
n = 5;

% parameter definition
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


% Initial guess definition
P = full(sprandsym(n,0.7,0.1,1));

options.maxiter = 1000;
options.ngtol = 1.0e-6;
options.eps2 = 1.0e-13;
 
 
options.a = 1;
options.b = 1;

        
I = eye(n);
        
% Riemannian gradient of f1 at P
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
[info] = algorithm3(P,options,theta);


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