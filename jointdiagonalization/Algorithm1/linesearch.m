function [alpha,flag] = linesearch(Y,V,MA,N,n,p,sigma,theta,alphatol)
flag = 0;
alpha = 1.0;  
phiY = meritfun(Y,MA,N,n,p);
while 1
    Q = qf(Y+alpha*V);
    phiQ = meritfun(Q,MA,N,n,p);
    test = phiQ-(1+2*sigma*theta*alpha)*phiY;
    if test > 0.0
        alpha = 0.5*alpha;
        if alpha < alphatol
            flag = 1;
            break;
        end
    else       
        break;
    end
end
end