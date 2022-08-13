function m1 = cmin(F)
%
% This function returns the minimum of a cell F
% It works toghether perf function only.
%
nalg = size(F,2);
m1 = min(F{1}{1});
for i = 2: nalg
  m2 =  min(F{i}{1});
  if m2 < m1
      m1 = m2;
  end
end

end