%
%   Data analysis
%  
%  Algorithm 1 x Algorithm 3 (JCAM)

A1 = readtable('echo1.dat');   % Algorithm 1
A3 = readtable('echo3.dat');   % Algorithm 3


h = semilogy(A1.iter,A1.gradnorm,'-or',...
             A3.iter,A3.gradnorm,'-sb');
legend('Algorithm 1 ','Algorithm 3','Location','Southwest')

xlabel('Iterate');
ylabel('Something...???');
set(h,'LineWidth',2); 

for i = 1 : size(A3,1)
    if i <= size(A1,1)
    fprintf("$%5d$ & $%15.10e$ & $%5d$  & NLS & $%15.10e $ &   $%5d$ \\\\ \n",A3.iter(i),A1.gradnorm(i),A1.evalf(i),A3.gradnorm(i),A3.evalf(i));
    else
    fprintf("$%5d$ &                    &          & NLS & $%15.10e $ &   $%5d$\\\\ \n",A3.iter(i),A3.gradnorm(i),A3.evalf(i)); 
    end
end