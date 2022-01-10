function echoiterinfo(etask,iecho,iter,nrGF,isVnonTS,alpha,dir)
if etask == 1
  fprintf(iecho,"  %5d   %25.20e %30.20e %20.10e %5s\n",iter,nrGF,isVnonTS,alpha,dir);
end
end