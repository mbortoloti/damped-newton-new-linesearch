function [sol] = algorithm3(P,options,theta)
    %
    %     Input description
    %
    %     options.maxiter ::: maximum of iterates
    %     options.ngtol   ::: gradient norm tolerance
    %     options. eps2   ::: tolerance for verification if solution
    %                         of Newton equation belongs to tangent space
    %     options.a       ::: constant a on f_1 definition
    %     options.b       ::: constant b on f_1 definition
    %     options.rgradf  ::: Riemannian gradient of f_1 definition
    %     options.metric  ::: Metric definition, as in (21)
    %     P               ::: initial guess
    %
    n = size(P,1);
    I = eye(n);   
    maxiter = options.maxiter;
    ngtol = options.ngtol;
    eps2 = options.eps2;
    ret = options.ret;
    rgradf = options.rgradf;
    metric = options.metric;
    
    sigma = 1.0e-3;
    a = options.a;
    b = options.b;
    
    % Gradient of merit function
    gphi = @(P) a*b*I-b^2*P^(-1); 
    
    % Merit function definition
    phi = @(P) 0.5* metric(rgradf(P),rgradf(P),P);
    
    sol.error = 0;
    
    dir = '   ';
    k = 0;
    alpha = NaN;
    
    
    % Loop for sequence generation
    fprintf('  Iterate       Gradient Norm       Step Length       Direction \n');
    tic;
    while k <= maxiter
        rgradfp = rgradf(P);
        ng = sqrt(metric(rgradfp,rgradfp,P));
        fprintf(' %5d  %23.12e %+18.12e %8s\n',k,ng,alpha,dir);
        if ng <= ngtol
            sol.time = toc;
            fprintf('\n ::: Gradient norm tolerance achieved\n');
            fprintf(' ::: Elapsed time  %12.8f\n',sol.time);
            break;
        end
        k = k +1;
        % Newton Equation (Sylvester equation)        
        sign = 1;
        try
        V = sylvester(P,P,2*(P^2-a/b*P^3));
        catch
        sign = 0;
        end
        
        
        gphip = gphi(P);

        if sign == 1
        % Verification if V belongs to tangent space
          if sqrt(metric(V-V',V-V',P)) <= eps2
            dir = 'new';
            if(snewtontest2(V,gphip,P,metric,theta) == 0)
                dir = 'gr2';
                V=-gphip;
            end
          else
            dir = 'grd';
            V = -gphip;
          end
        else
         V = -gphip;
        end
        
        % Armijo line search 
        alpha = armijo(P,V,gphip,sigma);
        
        % Iterate updating by exponential map
        P = ret(P,alpha*V);
    end
    if k<= maxiter
        sol.P = P;
    else
        sol.error = 1;
    end
    
    % Linesearch definition (Armijo type)
    function alpha = armijo(P,V,gphip,sigma)
        alpha = 1;
        phip = phi(P);
        while true
            phiq = phi(ret(P,alpha*V));
            arm = phiq - phip - alpha*sigma*metric(gphip,V,P);
            if arm <= 0
                break;
            else
                alpha = 0.5*alpha;
            end
        end
    end
    
    
end