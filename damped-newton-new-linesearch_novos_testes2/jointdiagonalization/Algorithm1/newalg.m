function [Y,info] = newalg(Y0,MA,N,n,p,options)
iecho = options.iecho;
etask = options.etask;
maxiter = options.maxiter;
sttol = options.sttol;
sigma = options.sigma;
theta = options.theta;
alphatol = options.alphatol;
tol = options.tol;
iter = 0;
Y = Y0;
isVnonTS = -1;
info.error = 0;
dir = '    ';
alpha = 0;
evalfa = NaN;
evalfn = NaN;

tstart = tic;
while 1
  %
  % Calculate Grad f(Y)
  %
  rGF = gradf(Y,MA,N,n,p);
  
  %
  % Calculate || Grad f(Y) ||
  %
  nrGF = norm(rGF,'fro');
   
  fprintf("  %5d   %25.20e %30.20e  %15.10e  %5s\n",iter,nrGF,isVnonTS,alpha,dir);
  echoiterinfo(etask,iecho,iter,nrGF,isVnonTS,alpha,dir,evalfa,evalfn);
  
  %
  %   evaluate function in Armijo linesearch
  %
  evalfa = 0;
  
  %
  %  evaluated function in new linesearch
  %
  evalfn = 0;
  
  if nrGF < tol
      info.time = toc(tstart);
      info.iter = iter;
      break;
  end
  iter = iter + 1;
  
  if iter > maxiter
      info.error = 1;
      info.time = -1;
      info.iter = iter;
      break;
  end
  
  % Compute Yo 
  Yo = comporth(Y);
  
  % Compute Zl 
  MZ = makeMZ(MA,Y,Y,N,n);
  
  % Compute Zlo
  MZo = makeMZ(MA,Y,Yo,N,n);
  
  % Compute Zlo2
  MZoo = makeMZ(MA,Yo,Yo,N,n);
  
  % Compute Right Hand Side of Newton Equation 
  nRHS = RHSnewton(MZ,MZo,N,n,p);
  
  % Compute Left Hand Side of Newton Equation
  nLHS = LHSnewton(MZ,MZo,MZoo,N,n,p);
  
  [Sn,nflag] = gmres(nLHS,-nRHS,size(nLHS,1),1.e-6,size(nLHS,1));
  if nflag > 0
      dir = 'grds';
      [Gmf] = gradmerit(Y,rGF,MA,N,n,p);
      V = -Gmf;
      [alpha,evalfa] = armijo(Y,V,Gmf,MA,N,n,p,sigma,evalfa);
  else
      B = invveck(Sn(1:p*(p-1)/2),p);
      C = invvec(Sn(p*(p-1)/2+1:p*(p-1)/2+p*(n-p)),n-p,p);
  
      V = Y*B+Yo*C;
      %
      % Tangent Space Test
      %
      isVnonTS = norm(Y'*V+V'*Y);
      if isVnonTS < sttol        
          [alpha,evalfn] = linesearch(Y,V,MA,N,n,p,sigma,theta,alphatol,evalfn);
          dir = 'new';
      else
          dir = 'grd';
          [Gmf] = gradmerit(Y,rGF,MA,N,n,p);
          V = -Gmf;
          [alpha,evalfa] = armijo(Y,V,Gmf,MA,N,n,p,sigma,evalfa);
      end
  end
  
  
  % Updating solution
  Y = qf(Y+alpha*V);
  
 
end

%info.iter = iter;
end