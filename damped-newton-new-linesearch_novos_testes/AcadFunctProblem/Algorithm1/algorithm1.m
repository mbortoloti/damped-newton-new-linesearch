function [P,sol] = algorithm1(P,options,theta)
    %
    %     Input description
    %
    %     options.maxiter ::: maximum of iterates
    %     options.ngtol   ::: gradient norm tolerance
    %     options. eps2   ::: tolerance for verification if solution
    %                         of Newton equation belongs to tangent space
    %     options.stpmin  ::: minimun step lenght
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
    stpmin = options.stpmin;
    
    % Evaluate number of the meriti function 
    evalf = 0;
    
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
        if ng < ngtol
            sol.time = toc;
            sol.iter = k;
            sol.evalf = evalf;
%             fprintf('\n ::: Gradient norm tolerance achieved\n');
%             fprintf(' ::: Elapsed time  %12.8f\n',sol.time);
%             fprintf(' ::: Elaluated function number  %5d\n',sol.evalf);
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
            [alpha,error,evalf] = newLS(P,V,sigma,stpmin,evalf);
            if error > 0
                V = -gphip;
                [alpha,~,evalf] = armijo(P,V,gphip,sigma,stpmin,evalf);
            end
          else
            dir = 'grd';
            V = -gphip;
            % Armijo line search 
            %alpha = armijo(P,V,gphip,sigma);
            [alpha,LSerror,evalf] = armijo(P,V,gphip,sigma,stpmin,evalf);
            if LSerror > 0
               sol.error = 3;
            end
          end
        else
            dir = 'grd';
            V = -gphip;
            % Armijo line search 
            %alpha = armijo(P,V,gphip,sigma);
            [alpha,LSerror,evalf] = armijo(P,V,gphip,sigma,stpmin,evalf);
            if LSerror > 0
                sol.error = 3;
            end
        end
        
        % Armijo line search 
        %alpha = armijo(P,V,gphip,sigma);
%         [alpha,LSerror,evalf] = armijo(P,V,gphip,sigma,stpmin,evalf);
%         if LSerror > 0
%             sol.error = 3;
%         end
%         
        % Iterate updating by exponential map
        P = ret(P,alpha*V);
    end
    if k > maxiter
        sol.error = 1;
    end
    
    % Linesearch definition (Armijo type)
    function [alpha,error,evalf] = armijo(P,V,gphip,sigma,stpmin,evalf)
        error = 0;
        alpha = 1;
%         phip = phi(P);
        [phip,evalf] = fevaluate(evalf,phi,P);
        GV = sigma * metric(gphip,V,P);
        while true
            %phiq = phi(ret(P,alpha*V));
            [phiq,evalf] = fevaluate(evalf,phi,ret(P,alpha*V));
            arm = phiq - phip - alpha * GV;
            if arm > 0
                alpha = alpha * 0.5;
                if alpha < stpmin
                    error = 1;
                    break;
                end
            else
                break;
            end
        end
    end
    
    % Linesearch definition (Armijo type)
    function [alpha,error,evalf] = newLS(P,V,sigma,stpmin,evalf)
        error = 0;
        alpha = 1.0;
%         phip = phi(P);
        [phip,evalf] = fevaluate(evalf,phi,P);
%         GV = sigma * metric(gphip,V,P);
        while true
            %phiq = phi(ret(P,alpha*V));
            [phiq,evalf] = fevaluate(evalf,phi,ret(P,alpha*V));
            arm = phiq - (1.0 + 2.0 * sigma * theta * alpha) * phip;
            if arm > 0
                alpha = alpha * 0.5;
                if alpha < stpmin
                    error = 1;
                    break;
                end
            else
                break;
            end
        end
    end
end