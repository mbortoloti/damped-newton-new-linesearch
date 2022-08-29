%
%   Data analysis
%  n = 100
%  First guess

TR = readtable('infoTR.txt');
BFGS = readtable('infoBFGS.txt');
ALGO1 = readtable('infoALGO1.txt');


h = semilogy(TR.iter,TR.gradnorm,'-or',...
         BFGS.iter,BFGS.gradnorm,'-sc',...
         ALGO1.iter,ALGO1.gradnorm,'-*b');
legend('Trust Region in [6]','BFGS in [6]','Algorithm 1','Location','Southwest')
xlabel('Iterate');
ylabel('|| grad f ||');
set(h,'LineWidth',2); 
set(gca,'FontSize',12);