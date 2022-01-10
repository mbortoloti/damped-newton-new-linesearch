function [Y,info] = damped(Y0,MA,N,n,p,options)
maxiter = options.maxiter;
etask = options.etask;
iecho = options.iecho;
% info.error = 1;
tol = options.tol;
sttol = options.sttol;
sigma = options.sigma;
iter = 0;
Y = Y0;
isVnonTS = 1;
dir = '    ';
alpha = 0;
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
   
  fprintf("  %5d   %25.20e %30.20e %20.10e %5s\n",iter,nrGF,isVnonTS,alpha,dir);
  
  echoiterinfo(etask,iecho,iter,nrGF,isVnonTS,alpha,dir);
  
  if nrGF < tol
      info.time = toc(tstart);
      info.iter = iter;
      info.error = 0;
      break;
  end
  iter = iter + 1;
  
  if iter > maxiter
      info.time=toc(tstart);;
      info.error = 1;
      info.iter = iter;
      break;
  end
  
  %
  % Grandient of Merit Function
  %
  Gmf = gradmerit(Y,rGF,MA,N,n,p);
  
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
  
  [Sn,sflag] = gmres(nLHS,-nRHS,size(nLHS,1),1.e-6,size(nLHS,1));
  if sflag > 0
      dir = 'grd';
      V = -Gmf;
  else
      B = invveck(Sn(1:p*(p-1)/2),p);
      C = invvec(Sn(p*(p-1)/2+1:p*(p-1)/2+p*(n-p)),n-p,p);
      V = Y*B+Yo*C;
      
      %
      % Tangent Space Test
      %
      isVnonTS = norm(Y'*V+V'*Y);
      if isVnonTS < sttol
          dir = ' new';
      else
          dir = ' grd';
          V = -Gmf;
      end
  end
      
  %
  %  Compute Lenght of Step 
  %
  alpha = armijo(Y,V,Gmf,MA,N,n,p,sigma);
  
  % Updating solution
  Y = qf(Y+alpha*V);
end


end