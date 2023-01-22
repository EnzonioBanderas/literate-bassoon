function func_rmTrialFromEvents(fp_d)
  
  fprintf('\n%s\n%s - Removing ''Trial'' events after eyetrack sync\n',...
    repmat('=',72,1),datestr(now,'yyyy-mm-dd hh:MM:ss'))


  fn_ptab = fullfile(fp_d,'participants.tsv');
  ptab = ter_readptab(fn_ptab);
  for i=1:size(ptab,1)
    pID = ptab.participant_id{i};
    fl = ter_listFiles(fullfile(fp_d,pID),[pID '*_events.tsv'],2);
    for j=1:numel(fl)
      fn_et = strrep(fl{j},'events.tsv','eyetrack.tsv');
      if exist(fn_et,'file')==2 || exist([fn_et '.gz'],'file')==2
        ev0 = ter_readBidsTsv(fl{j});
        ev = ev0(~ismember(lower(ev0.trial_type),'trial'),:);
        if not(isequaln(ev,ev0))
          fprintf('removing "Trial" events from : %s\n',fl{j});
          ter_writeBidsTsv(ev,fl{j});
        end
      end
    end
  end



end