function func_discardMp2rageFiles(fp_d,fp_dis)

  fl = ter_listFiles(fp_d,'*mp2rage*');
  
  fp_disc = fullfile(fp_dis,'mp2rage');
  
  
  for i=1:numel(fl)
    fn_old = fl{i};
    fn_new = strrep(fn_old,fp_d,fp_disc);
    fp_out = fileparts(fn_new);
    if ~isfolder(fp_out)
      mkdir(fp_out)
    end
    if exist(fn_new,'file')~=2
      movefile(fn_old,fn_new)
    end
  end