function rdtab = ter_rawDataInfo(fl_rawdata)

strPatBhdr = {...
  's' 'PatientName' '<ParamString."tPatientName">'
  };
for i=1:size(strPatBhdr,1)
  strPatBhdr{i,3} = regexprep(strPatBhdr{i,3},'\.','\\.');
  %strPatBhdr{i,3} = regexprep(strPatBhdr{i,3},'\<','\\<');
  %strPatBhdr{i,3} = regexprep(strPatBhdr{i,3},'\>','\\>');
  strPatBhdr{i,3} = strcat(strPatBhdr{i,3},'\ *\{[\ \n]*');
  %<ParamString."tPatientName">  { "Z7t2627^Fmrt Sfb1280a05study2"  }
  if strcmpi(strPatBhdr{i,1},'d')
    %strPatHdr{i,3} = strcat(strPatHdr{i,3},'>[\-0-9\.]+)');
  elseif strcmpi(strPatBhdr{i,1},'s')
    strPatBhdr{i,3} = strcat(strPatBhdr{i,3},'"(?<',strPatBhdr{i,2},...
      '>[\ \^a-zA-Z0-9\_\.\-]+)"');
  end
end

strPatHdr = {...
  's' 'ProtocolName'        'tProtocolName'
  's' 'idRefImage0'         'tReferenceImage0'  
  'd' 'baseRes'             'sKSpace.lBaseResolution'
  'd' 'nLinePE'             'sKSpace.lPhaseEncodingLines'
  'd' 'nSlices'             'sKSpace.lImagesPerSlab'
  'd' 'PosSag'              'sSliceArray.asSlice[0].sPosition.dSag'
  'd' 'PosCor'              'sSliceArray.asSlice[0].sPosition.dCor'
  'd' 'PosTra'              'sSliceArray.asSlice[0].sPosition.dTra'
  'd' 'PosAngleSag'         'sSliceArray.asSlice[0].sNormal.dSag'
  'd' 'PosAngleCor'         'sSliceArray.asSlice[0].sNormal.dCor'
  'd' 'PosAngleTra'         'sSliceArray.asSlice[0].sNormal.dTra'
  'd' 'thickness'           'sSliceArray.asSlice[0].dThickness'
  'd' 'FOVph'               'sSliceArray.asSlice[0].dPhaseFOV'
  'd' 'FOVro'               'sSliceArray.asSlice[0].dReadoutFOV'
  'd' 'TReff'               'sWiPMemBlock.alFree[9]'
  'd' 'kz'                  'sWiPMemBlock.alFree[18]'
  'd' 'ky'                  'sWiPMemBlock.alFree[19]'
  'd' 'TxRefVol'            'sTXSPEC.asNucleusInfo[0].flReferenceAmplitude'
  'd' 'AdjVol_sPosdSag'     'sAdjData.sAdjVolume.sPosition.dSag'
  'd' 'AdjVol_sPosdCor'     'sAdjData.sAdjVolume.sPosition.dCor'
  'd' 'AdjVol_sPosdTra'     'sAdjData.sAdjVolume.sPosition.dTra'
  'd' 'AdjVol_sNormdSag'    'sAdjData.sAdjVolume.sNormal.dSag'
  'd' 'AdjVol_sNormdCor'    'sAdjData.sAdjVolume.sNormal.dCor'
  'd' 'AdjVol_sNormdTra'    'sAdjData.sAdjVolume.sNormal.dTra'
  'd' 'AdjVol_dThickness'   'sAdjData.sAdjVolume.dThickness'
  'd' 'AdjVol_dPhaseFOV'    'sAdjData.sAdjVolume.dPhaseFOV'
  'd' 'AdjVol_dReadoutFOV'  'sAdjData.sAdjVolume.dReadoutFOV'
  'd' 'AdjVol_dInPlaneRot'  'sAdjData.sAdjVolume.dInPlaneRot'
  };

for i=1:size(strPatHdr,1)
  strPatHdr{i,3} = regexprep(strPatHdr{i,3},'\[','\\[');
  strPatHdr{i,3} = regexprep(strPatHdr{i,3},'\]','\\]');
  strPatHdr{i,3} = regexprep(strPatHdr{i,3},'\.','\\.');
  
  if strcmpi(strPatHdr{i,1},'d')
    %strPatHdr{i,3} = sprintf('%s\\ *= (?<%s>[\\-0-9\\.]+)',...
    %  strPatHdr{i,3},strPatHdr{i,2});
    strPatHdr{i,3} = strcat(strPatHdr{i,3},'\ *= (?<',strPatHdr{i,2});
    strPatHdr{i,3} = strcat(strPatHdr{i,3},'>[\-0-9\.]+)');
  elseif strcmpi(strPatHdr{i,1},'s')
    %strPatHdr{i,3} = sprintf('%s\\ *= (?<%s>[\\-0-9\\.]+)',...
    %  strPatHdr{i,3},strPatHdr{i,2});
    strPatHdr{i,3} = strcat(strPatHdr{i,3},'\ *= "(?<',strPatHdr{i,2});
    strPatHdr{i,3} = strcat(strPatHdr{i,3},'>[ \^a-zA-Z0-9\_\.\-\+]+)"');
  end
  strPatHdr{i,3} = regexprep(strPatHdr{i,3},'\\\\','\\');
end

for i=1:size(fl_rawdata,1)
  if size(fl_rawdata,2) ==1
    [~,text_hdr,text_bhdr] = ter_getSomeRawDataHeaderInfo(fl_rawdata{i});
  else
    text_hdr  = fileread(fl_rawdata{i,2});
    text_bhdr = fileread(fl_rawdata{i,1});
  end
  for j=1:size(strPatBhdr,1)
    tmp0 = regexp(text_bhdr,strPatBhdr{j,3},'names');
    fn = fieldnames(tmp0);
    if ~isempty(tmp0)
      if strcmpi(strPatBhdr{j,1},'d')
        tmp.(fn{1}) = str2double(tmp0(1).(fn{1}));
      else
        tmp.(fn{1}) = {tmp0(1).(fn{1})};
      end
    else
      if strcmpi(strPatBhdr{j,1},'d')
        tmp.(fn{1}) = 0;
      else
        tmp.(fn{1}) = {'NA'};
      end
    end  
  end
  for j=1:size(strPatHdr,1)
    tmp0 = regexp(text_hdr,strPatHdr{j,3},'names');
    fn = fieldnames(tmp0);
    if ~isempty(tmp0)
      if strcmpi(strPatHdr{j,1},'d')
        tmp.(fn{1}) = str2double(tmp0(1).(fn{1}));
      else
        tmp.(fn{1}) = {tmp0(1).(fn{1})};
      end
    else
      if strcmpi(strPatHdr{j,1},'d')
        tmp.(fn{1}) = 0;
      else
        tmp.(fn{1}) = {'NA'};
      end
    end  
  end
  tmp(1).ProtocolName{1} = strrep(tmp(1).ProtocolName{1},'+AF8-','_');
  date_tmp = strsplit(tmp.idRefImage0{1},'.');
  tmp.Studydate = str2double(date_tmp{end}(1:8));
  tmp = rmfield(tmp,'idRefImage0');
  if i==1
    rdtab = struct2table(tmp);
  else
    rdtab = [rdtab; struct2table(tmp)]; %#ok<AGROW>
  end
end
rdtab = rdtab(:,[1,end,2:end-1]);


%regexp(txt_bhdr{1},strPatBhdr{i,3},'names')
%tmp.ProtocolName = strrep(tmp.ProtocolName,'+AF8-','_');