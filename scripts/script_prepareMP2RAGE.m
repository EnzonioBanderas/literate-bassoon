fp_d   = '/media/diskEvaluation/Evaluation/sfb1280a05study7/dumpHereForSorting/rawdata';
fp_de  = '/media/diskEvaluation/Evaluation/sfb1280a05study7/dumpHereForSorting/derivatives';
fp_dis = '/media/diskEvaluation/Evaluation/sfb1280a05study7/dumpHereForSorting/discarded';
addpath(genpath('/media/diskEvaluation/Evaluation/sfb1280a05study7/scripts/packages/MP2RAGE-utils'))
addpath /media/diskEvaluation/Evaluation/sfb1280a05study7/scripts/packages/MP2RAGE-utils/core


fl_uni = essbids_listFiles(fullfile(fp_d,'**','anat','sub-*UNIT1*.nii.gz'));
fl_uni = fl_uni(~contains(fl_uni,'efaced'));
for i=1:numel(fl_uni)
  fp = fileparts(fl_uni{i});
  % find additional anat volumes to deface
  fl_add = essbids_listFiles(fullfile(fp,'sub-*.nii.gz'));
  fl_add = fl_add(not(ismember(fl_add,fl_uni{i})));
  %essbids_defaceAnatomicalVolume(fl_uni{i},fl_add,fp_dis);
end

nSlices = 256;
TRflash = 5.5e-3;  % echo spacing in seconds

txt0 = sprintf([...
  '# Required input parameters (full path must be provided)\n'...
  'INV1: FN_INV1\n'...
  'INV2: FN_INV2\n'...
  'UNI: FN_UNI\n'...
  '\n'...
  '# Optional input parameters (full path must be provided)\n'...
  'T1map:\n'...
  'sa2rageINV2:\n'...
  'sa2rageB1map:\n'...
  'tflB1map:\n'...
  '\n'...
  '# MP2RAGE parameters from protocol\n'...
  'B0: MY_B0 # Field Strength\n'...
  'TR: MY_TR # MP2RAGE TR\n'...
  'TRFLASH: MY_ECHOSPACING # GRE readout TR\n'...
  'FlipDegrees: MY_FAS # Flip angles\n'...
  'TIs: MY_INVTIMES # Inversion times\n'...
  'SlicesPerSlab: MY_SLICES\n'...
  'PartialFourierInSlice: MY_PF\n'...
  '\n'...
  '# Background removal options (1 = yes, 0 = no)\n'...
  'DenoiseUNI: 1\n'...
  'DenoiseT1map: 1\n'...
  'DenoiseWeight: 10\n'...
  '\n'...
  '# Calculate additional maps (1 = yes, 0 = no)\n'...
  'CalculateT1map: 1\n'...
  'CalculateR1map: 0\n'...
  'CalculateM0map: 0\n'...
  '\n'...
  '# B1+ correction (1 = yes, 0 = no), coregister requires SPM\n'...
  'B1correct: 1\n'...
  'Coregister: 0\n'...
  '\n'...
  '# If B1+ correcting with Sa2RAGE, not used here\n'...
  'sa2rageTR: 2.4\n'...
  'sa2rageTRFLASH: 2.75e-3\n'...
  'sa2rageFlipDegrees: 4 10\n'...
  'sa2rageTIs: 47e-3 1800e-3\n'...
  'sa2rageBaseResolution: 128\n'...
  'sa2ragePartialFourierInPE: 0.75\n'... 
  'sa2rageiPATPhaseEncode: 3\n'...
  'sa2rageRefLines: 18\n'...
  'sa2rageAverageT1: 1.5\n']);
% activate spm if not yet started up
try 
  spm('version');
catch
  startspm 12;
  spm quit;
end
fl_uni = [
  essbids_listFiles(fullfile(fp_d,'**','sub-*UNIT1*.nii.gz'));
  essbids_listFiles(fullfile(fp_dis,'**','sub-*UNIT1*.nii.gz'))];
fl_uni = fl_uni(~contains(fl_uni,'efaced'));
for i=1:numel(fl_uni)
  fn = fullfile(fl_uni{i});
  bf = essbids_parseLabel(fn);
  fp_out = fullfile(fp_de,'MP2RAGE_utils',['sub-' bf.sub],['ses-' bf.ses]);
  if not(isfolder(fp_out))
    mkdir(fp_out);
  end
  fl_clean = dir(fullfile(fp_out,'*clean*.nii*'));
  if numel(fl_clean)>0
    fprintf('Already denoised and T1map calculated, skipping %s\n',fp_out);
    continue;
  end
  bf1 = bf;
  bf1.suffix = 'MP2RAGE';
  bf1.inv = '1';
  fn0_inv1 = fullfile(bf1.fpath,essbids_buildFileName(bf1));
  bf1.inv = '2';
  fn0_inv2 = fullfile(bf1.fpath,essbids_buildFileName(bf1));
  bf1 = essbids_parseLabel(fn0_inv1);
  bf2 = essbids_parseLabel(fn0_inv2);
  jinfo1 = jsondecode(fileread(essbids_findJsonSidecar(fn0_inv1)));
  jinfo2 = jsondecode(fileread(essbids_findJsonSidecar(fn0_inv2)));
  fn1_uni = fullfile(fp_out,[bf.fname '.nii']);
  fn1_inv1 = fullfile(fp_out,[bf1.fname '.nii']);
  fn1_inv2 = fullfile(fp_out,[bf2.fname '.nii']);
  if exist(fn1_uni,'file')~=2
    if contains(lower(bf.extension),'.gz')
      gunzip(fn,fp_out);
      gunzip(fn0_inv1,fp_out);
      gunzip(fn0_inv2,fp_out);
    else
      copyfile(fn,fp_out);
      copyfile(fn0_inv1,fp_out);
      copyfile(fn0_inv2,fp_out);
    end
  end
  fn_config = fullfile(fp_out,['sub-' bf.sub '_config.yml']);
  %delete(fn_config)
  if exist(fn_config,'file')~=2
    txt = txt0;
    txt = strrep(txt,'MY_B0',sprintf('%.1f',jinfo1.MagneticFieldStrength));
    txt = strrep(txt,'FN_UNI', fn1_uni);
    txt = strrep(txt,'FN_INV1',fn1_inv1);
    txt = strrep(txt,'FN_INV2',fn1_inv2);
    txt = strrep(txt,'MY_TR',sprintf('%.3f',jinfo1.RepetitionTime));
    txt = strrep(txt,'MY_FAS',sprintf('%d %d',...
      jinfo1.FlipAngle,jinfo2.FlipAngle));
    txt = strrep(txt,'MY_INVTIMES',sprintf('%.3f %.3f',...
      jinfo1.InversionTime,jinfo2.InversionTime));
    txt = strrep(txt,'MY_PF',sprintf('%d',jinfo1.PartialFourier));
    txt = strrep(txt,'MY_SLICES',sprintf('%d',nSlices));
    txt = strrep(txt,'MY_ECHOSPACING',sprintf('%d',TRflash));
    fid = fopen(fn_config,'w');
    fprintf(fid,txt);
    fclose(fid);
  end
  mp2rage_main(fn_config);
  eval(sprintf('!pigz -9 %s',fullfile(fp_out,'*.nii')));
end


