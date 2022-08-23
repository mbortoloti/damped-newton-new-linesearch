function [fX,eval] = fevaluate(eval,f,X)
    eval = eval + 1;
    fX = f(X);
end