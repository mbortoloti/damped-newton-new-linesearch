%
%  Driver program for Joint Diagonalization
%  Newton Method

nG = 10; % Number of Initial Guesses

%Echo task 
% etask = 0 --> do not print echo info
% etask = 1 --> do print echo info   
etask = 1;

% Files setting
iecho = fopen('echo.dat','w');   % Echo file
itime = fopen('natime.dat','w');   % Performance profile 


%N = 10;
% dim_n = [5,10, 10, 20, 20, 25, 25];
% dim_p = [3, 3,  5,  5, 10, 10, 15];
N = 5; 
dim_n = [ 10,10,10,10];
dim_p = [  2, 3, 4, 5];
rho = [0.0001,0.001,0.01,0.1,1,10,100];
pscat = zeros(size(rho));
% dim_n = [10,50,100,300];
% dim_p = [2, 10, 20, 50];



% address = pwd;
address = '/home/marcio/GDrive/TelesMarcio/MatlabCodes/Jointdiagonalization';

formatSpec = '%18.12f';
adguess = strcat(address,'/guesses/');

%
% Newton options card
%
% options.tol = 1.e-12;
% options.maxiter = 5000;
options.iecho = iecho;
  options.maxiter    = 5000;
    options.tol      = 1.e-06;
    options.sttol    = 1.e-12;
    options.sigma    = 1.e-4;
    options.theta    = 0.9;
    options.alphatol = 1.e-4;   
    options.etask = etask;
    
for k = 1:size(dim_n,2)
  n = dim_n(k);
  p = dim_p(k);
  
  echodim(iecho,n,p,etask);
  
%   smname = 'mA_';
%   smname = strcat(smname,int2str(p),'_',int2str(n),'.dat');
%   sfile = strcat(adguess,smname);
%   A = getinput(sfile,formatSpec,p,n);
%   
%   fileecho(iecho,smname,etask);
  
  smname = 'mMA_';
  smname = strcat(smname,int2str(n),'_',int2str(n),'.dat');
  sfile = strcat(adguess,smname);
  MA = getinput(sfile,formatSpec,n,N*n);
  
  fileecho(iecho,smname,etask);
  
%   smname = 'mMC_';
%   smname = strcat(smname,int2str(p),'_',int2str(p),'.dat');
%   sfile = strcat(adguess,smname);
%   MC = getinput(sfile,formatSpec,p,m*p);
%   
%   fileecho(iecho,smname,etask);
  for kk = 1:size(rho,2)
  for i = 1:nG  % Loop for initial guess test
    %
    % Matrix A Reading Procedure
    %
    smname = 'guess_';
    smname = strcat(smname,int2str(n),'_',int2str(p),'_',int2str(kk),'_',int2str(i),'.dat');
    sfile = strcat(adguess,smname);
    Y0 = getinput(sfile,formatSpec,n,p);
 
    fileecho(iecho,smname,etask);
   
    [Y,info] = newalg(Y0,MA,N,n,p,options);
%     [Y,info] = newton(Y0,MA,N,n,p,options);
    
    
    infoecho(iecho,info,etask);
    [pscat] = mvtime(itime,info,kk,pscat);
  end
  end
  
end


fclose(iecho);
fclose(itime);



