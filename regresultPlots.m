cdir = loaddirfun;
mimsInfo = MIMSsummary;
anvec = [84,46,6,15,17,45,47,37,36]; %[3,7,44];
mimslabels = cellfun(@(x) x(10:end),mimsInfo.dsetInfo(anvec,1),'uniformoutput',0);
regData = cellstr(ls(cdir.saMIMSbrain));
for i = 1:7
  associatedfiles = find(cellfun(@(x) ~isempty(strfind(x,strtok(mimslabels{i},'_'))),regData));
  labelSubset = regData(associatedfiles);
  load([cdir.saMIMSbrain mimslabels{i} '.mat']); %loads dp, chanInd, I0, Jmat, usedat
  Ldat = load([cdir.saMIMSbrain mimslabels{i} '_Left.mat']); %loads Ldat.{'Mmimstsu', 'mu_tst'}
  Rdat = load([cdir.saMIMSbrain mimslabels{i} '_Right.mat']); %loads Rdat.{'Mmimstsu', 'mu_tst'}
  %rpM = regionprops(Ldat.Mmimstsu > 0,'boundingbox');
  
  Lreg = 
end