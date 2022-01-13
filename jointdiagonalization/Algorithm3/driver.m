clear all;
%
%  Driver for Newton 
%  Joint Diagonalization
%
N = 10;
n = 30;
p = 10;

%Echo task 
% etask = 0 --> do not print echo info
% etask = 1 --> print echo info   
etask = 0;

P = qf(randn(n,n));

d = sort(rand(n,1),'descend');
ML = diag(d);
MA = P*ML*P';
for l = 2:N
  d = sort(rand(n,1),'descend');
  MA = [MA,P*diag(d)*P'];
%   ML = [ML,diag(d)];
end

% Exact solution
Yex = P*eye(n,p);

% Initial guess
rho = 10;
Y0 = qf(Yex+rho*rand(n,p));

options.maxiter = 10000;
options.tol = 1.e-12;
options.sttol = 1.e-15;
options.sigma = 1.e-3;
options.etask = etask;
options.iecho = fopen("echo.dat","w");   %Echo file
[Y,info] = damped(Y0,MA,N,n,p,options);
