fn_newDat = '/media/diskEvaluation/Evaluation/sfb1280a05study7/sourcedata/mri_rawdata/meas_MID00429_FID14417_MP2RAGE_anat_ses_1_state_wip.dat';
%'/media/diskEvaluation/Evaluation/sfb1280a05study7/sourcedata/mri_rawdata/sub-Z7T7221_meas_MID00160_FID15731_MP2RAGE_anat_ses_1.dat';
fn_oldDat = '/media/diskEvaluation/Evaluation/sfb1280a05study2b/sourcedata/sub-EL07NZ11/raw_mp2rage/sub-EL07NZ11_meas_MID96_anat_mp2rage_FID60320.dat';


tic
[rdhdr1,text_hdr1,text_bhdr1] = ter_getSomeRawDataHeaderInfo(fn_oldDat);
toc
tic
[rdhdr2,text_hdr2,text_bhdr2] = ter_getSomeRawDataHeaderInfo(fn_newDat);
toc


fid = fopen(fn_newDat,'r');
tic
t         = fgets(fid);
toc
isheader = 0;
t2 = '';
while isheader<3
  t         = fgets(fid);
  if contains(t,'ASCCONV BEGIN')
    %contains(t,'ASCCONV BEGIN')
    %contains(t,'### ASCCONV BEGIN ###')
    isheader = 1;
    %i = strfind(t,'ASCCONV BEGIN');
    %i = strfind(t,'### ASCCONV BEGIN ###');
    %t = t(i:end);
  elseif contains(t,'ASCCONV END')
    isheader = 3;
  end
  if isheader>0
    disp(t);
  else
    t2 = [t2,t];
  end
end
fclose(fid);

disp(t2)
