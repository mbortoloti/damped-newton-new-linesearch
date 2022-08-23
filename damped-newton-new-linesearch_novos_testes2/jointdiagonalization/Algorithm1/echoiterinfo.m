function echoiterinfo(etask,iecho,iter,nrGF,isVnonTS,alpha,dir,evalfa,evalfn)
if etask == 1
   if strcmp('new',dir)
       fprintf(iecho,"%5d,%25.20e,%30.20e,%20.10e,%5s,%5d\n",iter,nrGF,isVnonTS,alpha,dir,evalfn);
   else
       fprintf(iecho,"%5d,%25.20e,%30.20e,%20.10e,%5s,%5d\n",iter,nrGF,isVnonTS,alpha,dir,evalfa);
   end
end
end