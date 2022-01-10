function [pscat] = mvtime(itime,info,kk,pscat)
  if info.error > 0
      fprintf(itime,'%25s\n','INF');
  else
      fprintf(itime,'%25.15f\n',info.time);
      pscat(kk)=pscat(kk)+1;
  end

end