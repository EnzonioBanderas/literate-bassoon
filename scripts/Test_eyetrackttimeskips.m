
fn='sub-Z7T7104_ses-1_run-1_datetime-20211203T144912_eyetrack.txt';
fn='sub-Z7T7104_ses-1_run-2_datetime-20211203T145651_eyetrack.txt';
fp='/media/diskEvaluation/Evaluation/sfb1280a05study7/pilotdata/pilot_2021-12-03_Z7T7104/ET';

mytxt=fileread(fullfile(fp,fn));

tmp=strsplit(mytxt,'\n');
ind = cellfun(@(x) regexp(x,'^10\t')==1 ,tmp ,'UniformOutput',0);
ind = not(cellfun(@isempty,ind));
tmp=tmp(ind);

fn_tmp = '/media/diskEvaluation/Evaluation/sfb1280a05study7/tmp.txt';
fid = fopen(fn_tmp,'w');
cellfun(@(x) fprintf(fid,'%s\n',x),tmp);
fclose(fid);

%t1=
t2=readtable(fn_tmp,'delim','tab');
delete(fn_tmp)


tmp=strsplit(mytxt,'\n');
ind = cellfun(@(x) regexp(x,'^12\t')==1 ,tmp ,'UniformOutput',0);
ind = not(cellfun(@isempty,ind));
tmp=tmp(ind);
%ev1=tmp;
ev2=tmp;

