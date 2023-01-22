function func_prepareCardio(fp_d,fp_de)
  
  fprintf('\n%s\n%s - Preparing data for cardio evaluation\n',...
    repmat('=',72,1),datestr(now,'yyyy-mm-dd hh:MM:ss'))
  
  fp_out = fullfile(fp_de,'cardioEval');
  ter_bidsPhysio2cardioEval(fp_d,fp_out);
  
end