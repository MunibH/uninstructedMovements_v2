function meta = loadMAH21_MCStim(meta,datapth)


meta(end+1).datapth = datapth;
meta(end).anm = 'MAH21';
meta(end).date = '2023-09-07';
meta(end).datafn = findDataFn(meta(end));
meta(end).datapth = fullfile(meta(end).datapth,'DataObjects',meta(end).anm,meta(end).datafn);
meta(end).stim = 'bilateral';
meta(end).stimLoc = 'Bi_MC';
meta(end).stimPow = 1; % mw/mm^2 (allpow=8)
meta(end).stimEpoch = 'delay'; 

% meta(end+1).datapth = datapth;
% meta(end).anm = 'MAH21';
% meta(end).date = '2023-09-12';
% meta(end).datafn = findDataFn(meta(end));
% meta(end).datapth = fullfile(meta(end).datapth,'DataObjects',meta(end).anm,meta(end).datafn);
% meta(end).stim = 'bilateral';
% meta(end).stimLoc = 'Bi_MC';
% meta(end).stimPow = 1; % mw/mm^2 (allpow=8)
% meta(end).stimEpoch = 'delay'; 
% 
% meta(end+1).datapth = datapth;
% meta(end).anm = 'MAH21';
% meta(end).date = '2023-09-13';
% meta(end).datafn = findDataFn(meta(end));
% meta(end).datapth = fullfile(meta(end).datapth,'DataObjects',meta(end).anm,meta(end).datafn);
% meta(end).stim = 'bilateral';
% meta(end).stimLoc = 'Bi_MC';
% meta(end).stimPow = 1; % mw/mm^2 (allpow=8)
% meta(end).stimEpoch = 'delay'; 


end

