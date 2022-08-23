%
% This code helps in the analyses of Damped method for 
%              non-conservative problems 
%
clear all;
clc;

nig = 1; % Number of Initial Guesses

options.maxiter = 5000;
options.ngtol = 1.0e-6;

% Dimension definition
dim = [500];
rng(12345,'twister');

itime = fopen("nlstime.dat","w");
ieft  = fopen("nlsift.dat","w");
iecho = fopen("echo.log","w");
options.iecho = iecho;

% Define seed for rand
rng(5000);

for n = dim

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

  for i = 1 : nig
    p = rand(n,1);
    p = p/norm(p);
    
    
    fprintf(iecho," Dimension  %5d  Initial guess  %5d\n",n,i);
    % Solver call
    [info] = nonconservativedamped(p,options);
     if info.error > 0
        fprintf(itime,'%25s\n','inf');
        fprintf(ieft ,'%25s\n','inf');
    else
        fprintf(itime,'%25.15f\n',info.time);
        fprintf(ieft ,   '%25d\n',info.eft);
    end    
  end
end
fclose(itime);
fclose(ieft);