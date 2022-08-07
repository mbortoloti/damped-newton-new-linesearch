%
%   logdet problem
%





clear all;
clc;

%
% Files for performance analysis
%
itime = fopen("logdettimeBFGS.dat","w");
% iiter = fopen("logdetiterates.dat","w");
% ieval = fopen("logdetfeval.dat","w");


%
% Static random numbers (only for tests)
%
rng(12345,'twister');

% Dimension definition 
dimensions = [100,200,300,400,500];
% dimensions = [5];
% Number of initial guesses for each dimension
nig = 10;

for nn = 1 :  size(dimensions,2)
n = dimensions(nn);
fprintf("n = %5d\n",n);

%manifold setting
manifold = sympositivedefinitefactory(n);
problem.M = manifold;


%
%  objective function
%
a = 1.0;
b = 1.0;
f = @(X) a * log(det(X)) + b * trace(X^-1);
problem.cost = f;

problem.egrad = @(X)(a * X^-1 - b * X^-2);

% checkgradient(problem);

for q = 1 : nig

% Initial guess definition
%P = full(sprandsym(n,0.7,0.1,1));
P0 = rand(n,n);
P0 = 0.5 * (P0 + P0') + n * eye(n);


 options.maxiter =5000;
 options.minstep = 1.e-10;
 options.tolgradnorm = 1.0e-6;

 
% Solver call 
%[P,info] = algorithm1(P,options,theta);
[P,cost,info,~] = rlbfgs(problem,P0,options);
 iterations = [info.iter]; 
 etime = [info.time];
 if iterations(end) >= options.maxiter
     fprintf(itime,"%20s\n","INF");
%     fprintf(iiter,"%10s\n","INF");
%     fprintf(ieval,"%10s\n","INF");
 else
     fprintf(itime,"%20.15e\n",etime(end));
%     fprintf(iiter,"%10d\n",info.iter);
%     fprintf(ieval,"%10d\n",info.evalf);
 end

% end of initial guess loop
end

%end of dimensions loop
end

fclose(itime);
% fclose(iiter);
% fclose(ieval);
