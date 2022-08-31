clear all;
%
%  Driver for Newton 
%  Joint Diagonalization
%
N = 5;
n = 20;
p = 9;

%Echo task 
% etask = 0 --> do not print echo info
% etask = 1 --> print echo info   
etask = 1;

%Static random generation
rng(12345,'twister');

% Setting MA matrices. Note that MA = [A_1 A_2 ... A_N]
P = qf(randn(n,n));
d = sort(rand(n,1),'descend');
ML = diag(d);
MA = P*ML*P';
for l = 2:N
  d = sort(rand(n,1),'descend');
  MA = [MA,P*diag(d)*P'];
end

% Exact solution
Yex = P*eye(n,p);

% Initial guess
rho = 0.1;
Y0 = qf(Yex+rho*rand(n,p));

options.maxiter = 10000;                 % Maximum of iterations
options.tol = 1.e-12;                    % Tolerance for norm of gradient function
options.sttol = 1.e-12;                  % Tolerance for tangent space test for solution of Newton's equation
options.sigma = 1.e-3;                   % Parameter for linesearch
options.etask = etask;                   % print echo file 
options.iecho = fopen("echo.dat","w");   %Echo file
[Y,info] = damped(Y0,MA,N,n,p,options);
