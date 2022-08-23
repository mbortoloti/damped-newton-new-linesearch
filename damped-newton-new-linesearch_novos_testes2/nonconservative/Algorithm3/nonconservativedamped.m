function   [info] = nonconservativedamped(p,options)
    % 
    %     "Damped Newton Method on Riemannian Manifolds" by
    %            Marcio Ant�nio de A. Bortoloti
    %            Teles A. Fernandes
    %            Orizon Ferreira
    %            Yuan J. Yun
    %
    %    This function finds the critical points of a 
    %    non-conservative field 
    % 
    %
    %     Input description
    %
    %     options.Q       ::: A skew-symmetric matrix
    %     options.pstar   ::: Critical point definition
    %     options.X       ::: Non-conservative function definition
    %     options.dX      ::: Riemannian derivative of X definition
    %     options.phi     ::: Merit function definition
    %     options.retr    ::: Exponential map definition
    %     options.maxiter ::: maximum of iterates
    %     options.ngtol   ::: gradient norm tolerance
    %     p               ::: initial guess
    %
    theta = 0.9999;
    iecho = options.iecho;
    sigma = 1.e-3;
    eps2 = 1.0E-08;
    eps1 = options.ngtol;
    R = options.retr;
    phi = options.phi;
    Q = options.Q;
    n = size(Q,1);
    I = eye(n);
    pstar = options.pstar;
    X = options.X;
    dX = options.dX;
    kMAX = options.maxiter;
    DIR = '   ';
    alpha = -1;
    k = 0;
    
    ng = norm(X(p));
    ef = 0;
    info.eft = 0;
    info.error = 0;
    fprintf('  Iterate     Gradient norm        Step length      Direction    F. eval.\n');
    fprintf(iecho,'  Iterate     Gradient norm        Step length      Direction    F. eval.\n');
    % Iterate loop
    tic;
    while true
        if k > kMAX
            fprintf(' ::: Maximum iterate number achieved.\n');
            fprintf(iecho,' ::: Maximum iterate number achieved.\n');
            info.iter = k;
            info.time = toc;
            info.p = p;
            info.ngradf = ng;
            info.error = 1001;
            break;
        end
            
        fprintf('%6d %23.12e %+18.12e %7s  %9d\n',k,ng,alpha,DIR,ef);
        fprintf(iecho,'%6d,%23.12e,%+18.12e,%7s,%9d\n',k,ng,alpha,DIR,ef);
        if ng <= eps1 
            info.time = toc;
            info.p = p;
            info.iter = k;
            info.ngradf = ng;
            fprintf('\n ::: Critical point found.\n');
            fprintf(iecho,'\n ::: Critical point found.\n');
            fprintf(' ::: Evaluate Time %12.5f\n',info.time);
            fprintf(iecho,' ::: Evaluate Time %12.5f\n',info.time);
            break;
        end
        ef = 0;
        k = k + 1;
        K = (I-p*p')*dX(p);
        Ka = [K;p'];
        b = -X(p);
        ba = [b;0];    
        
        % Gradient of merit function
        Gphip = (I-p*p')*dX(p)'*X(p);
              
        % Verification if Newton equation has a solution
        rankK = rank(Ka,eps2);
        rankb = rank([Ka,ba],eps2);
        if rankK == rankb
            DIR = 'new';
            % Newton equation solution
            v = Ka\ba;
            if Gphip'*v + theta * norm(Gphip)*norm(v) > 0
                DIR = 'grd';
                v = -Gphip;
            end
        else     
            DIR = 'grd';
            v = -Gphip;
        end
        
        % Armijo line search
        [alpha,ef] = armijo(p,v,Gphip,ef);
        if alpha < 1.e-10 
                info.error = 2;
                break;
        end
        info.eft = info.eft + ef;
        
        % Sequence updating
        p = R(p,alpha*v);
        
        ng = norm(X(p));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function [alpha,ef] = armijo(p,v,Gphip,ef)
        alpha = 1;
        phip = phi(p);
        ef = ef + 1;
        while true            
            q = R(p,alpha*v);
            arm = phi(q)-phip-sigma*alpha*Gphip'*v;
            ef = ef + 1;
            if arm <= 0
                break;
            else
                alpha= 0.5*alpha;
            end
        end
    end
end