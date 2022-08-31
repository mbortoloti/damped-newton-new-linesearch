%
%  A modified Armijo for damping a Riemannian Newton-type method
%  by Bortoloti, M. A. de A. and Fernandes, T. A.
%
%  Code for finding a zero of a skew-symmetric vector field 
%              
%
clear all;
clc;

% nig = 15; % Number of Initial Guesses

options.maxiter = 2000;
options.ngtol = 1.0e-6;




% Dimension setting
n = 5;

% A random skew-symmetric matrix generation
Q = rand(n,n);
Q = 0.5*(Q-Q');
options.Q = Q;
    
I = eye(n);
    
% Critical point definition for tests
pstar = ones(n,1);
pstar = pstar/norm(pstar);
options.pstar = pstar;

% Non-conservative field definition
X = @(p) (I-p*p')*Q*(p-pstar);
options.X = X;

% Field derivative definition
dX = @(p) (I+p*p')*Q-(p*(Q*(p-pstar))'+p'*Q*(p-pstar)*I);
options.dX = dX;

% Merit function definition
phi = @(x) 0.5*norm(X(x))^2;
options.phi = phi;

% Retraction Definition
% R = @(p,v)  p*cos(norm(v))+sin(norm(v))*v/norm(v);
R = @(p,v) (p+v)/norm(p+v);
options.retr = R;



% Initial guess
    p = rand(n,1);
    p = p/norm(p);
    
    

    % Solver call
    [info] = algorithm1(p,options);
