function infoecho(iecho,info,etask)
  if etask == 1
  fprintf(iecho,' ::: Solution Info\n');
  fprintf(iecho,' ::: Error code .................. %10d\n',info.error);
  fprintf(iecho,' ::: Elapsed time (s) ............ %10.5f\n',info.time);
  fprintf(iecho,' ::: Number of iterations ........ %10d\n',info.iter);
  fprintf(iecho,'\n');
  end
end

