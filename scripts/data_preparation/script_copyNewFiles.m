script_init_study6;

fp_sr = fullfile(fp0,'sourcedata_raw');
fl = ter_listFiles(fp_sr,'*');
fl = fl(not(isfolder(fl)));
t1 = struct2table(cellfun(@dir,fl));
t1 = sortrows(t1,'datenum');

t2 = t1(t1.datenum > datenum('2021-10-18 00:00:00'),:);
%t2 = t1(t1.datenum > datenum('2021-09-08 12:00:00'),:);
%t2 = t1(t1.datenum > datenum('2021-08-29 00:00:00'),:);
%t2 = t1(t1.datenum > datenum('2021-08-17 13:00:00'),:);
%t2 = t1(t1.datenum > datenum('2021-08-09 00:00:00'),:);
%t2 = t1(t1.datenum > datenum('2021-07-25 00:00:00'),:);
%t2 = t2(t2.datenum < datenum('2021-07-04 00:00:00'),:);
%t2=t1(contains(t1.name,'04DA18_'),:);
disp(t2);
fprintf(['%d files found\n'...
  'Press any key to continue or Crtl+c to abort\n'],size(t2,1));
pause
fprintf('Copying %d files ...\n',size(t2,1));

for i=1:size(t2,1)
  fn_old = fullfile(t2.folder{i},t2.name{i});
  fn_new = fullfile(fp_us,t2.name{i});
  if regexp(t2.name{i},'EKG_[a-zA-Z0-9]*\.pdf')>0
    fn_new = fullfile(fp_us, ['sub-' ...
      t2.name{i}(5:(strfind(t2.name{i},'.pdf')-1)) '_ecgscan.pdf']);
  end
  if exist(fn_new,'file')~=2
    copyfile(fn_old,fn_new)
  end
end

% remove spaces in file names
fl_withSpaces = ter_listFiles(fp_us,'* *');
for i=1:numel(fl_withSpaces)
  fn_old = fl_withSpaces{i};
  [fp,fn,fe] = fileparts(fn_old);
  fn_new = fullfile(fp,[strrep(fn,' ','') fe]);
  if not(exist(fn_new,'file')==2)
    movefile(fn_old,fn_new)
  end
end
fprintf('done.\n');
