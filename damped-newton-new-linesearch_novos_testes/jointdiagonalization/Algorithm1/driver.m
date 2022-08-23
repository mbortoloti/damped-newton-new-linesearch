clear all;
%
%  Driver for Algorithm 1 
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

% Exact solution for numerical tests
Yex = P*eye(n,p);

% Initial guess
rho = 0.1;
Y0 = qf(Yex+rho*rand(n,p));


% Parameters for code
options.maxiter = 10000;                      % Maximum number of iterations
options.tol = 1.e-12;                         % Tolerance for gradient norm of the function f
options.sttol = 1.e-12;                       % Tolerance for tangent space test of Newton's equation solution  
options.sigma = 1.e-3;                        % Parameter for linesearch
options.theta = 1e-3;                         % Parameter for linesearch
options.alphatol = 1.e-2;                     % Tolerance for alpha choice
options.iecho = fopen("echo.dat","w");        % Echo file
options.etask = etask;                        % Echo file print option
[Y,info] = newalg(Y0,MA,N,n,p,options);