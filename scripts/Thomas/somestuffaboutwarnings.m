




%warnstate = warning('query','MATLAB:datetime:FormatConflict_mM');
%warning('off','MATLAB:datetime:FormatConflict_mM')
% my code 
%warning(warnstate.state,'MATLAB:datetime:FormatConflict_mM')

bWarned = true;
warncode=5;
if bWarned
  warning('myprog:label:beWarned','Be warned : %05.2f',warncode)
end
[~,warnlabel] = lastwarn();
warnstate = warning('query',warnlabel);
warning('off',warnlabel)
warning('myprog:label:beWarned','Be warned')
warning(warnstate.state,warnlabel)
warning('myprog:label:beWarned','Be warned')

error('Ohno')