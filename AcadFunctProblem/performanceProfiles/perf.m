function [x,y,pf]=perf(varargin)
%
%   Input arguments (arq1,...,arqn,label1,...,labeln)
%
  nalg = (nargin)/2;
  for i = 1 : nalg
%       varargin{1,i}
    fid = fopen(varargin{1,i},'r');
    F{i} = textscan(fid,'%20.15f');
    fclose(fid);
    clear fid;
  end 
  npb = size(F{1}{1},1);
  for k = 1 : npb
      for kk = 1 : nalg
          P(kk) = F{kk}{1}(k);
      end
      m = min(P);
      for kk = 1 : nalg
          F{kk}{1}(k) = F{kk}{1}(k)/m;
      end
  end
%
% This procedure finds the cases could not reach the solution
% inside of step limit
%
%   for k = 1 : nalg   % Algorithm setting
%       ifail = 0;
%       for n = 1 : npb      % Problem setting
%         if F{k}{1}(n) == Inf
%             ifail = ifail + 1;
%         end
%       end
%       fail(k) = ifail;
%     end
 %
 % The pf matrix have the performance profile data
 %
  t = 0;
  pf = zeros(nalg,npb*nalg);
  %x = zeros(npb*nalg);
  while true
    m = cmin(F);
    if m == Inf 
        break;
    end
    t = t + 1;
    x(t) = m;
    for k = 1 : nalg   % Algorithm setting  
      for n = 1 : npb      % Problem setting
        if F{k}{1}(n) == m
          pf(k,t) = pf(k,t) + 1;
          F{k}{1}(n) = Inf;
        end     
      end
    end
  end
  
% This code sectin acumulates data in pf matrix in order to construct the
% performance profile plot
 
 y = zeros(nalg,size(x,2));
 
 for k = 1 : nalg   % Algorithm setting  
     y(k,1) = pf(k,1);
 end
 
%  y(:,1) = pf(:,1);
 
 for k = 1 : nalg   % Algorithm setting  
      for n = 2 : size(x,2)      % Problem setting
         y(k,n) = pf(k,n)+y(k,n-1); 
      end
 end
 
 y = y/npb;
color = ['c','b','g','r','y','m','k'];
mark  = ['o','x','*','+'];
%figure
% graphics_toolkit('fltk');
x = log10(x);
for k = 1:nalg
%    stairs(x,y(k,:),color(k),'LineWidth',2,'Marker',mark(k));  
 stairs(x,y(k,:),color(k),'LineWidth',2);     
    hold on;    
end
hold off;
for k = nalg+1 : nargin
    leg{k-nalg} = varargin{1,k};
end
%axis([0 0.72 0.2 1]);
legend(leg,'Location','SOUTHEAST');
saveas(gcf,'newtonxdampedall.png')
end
