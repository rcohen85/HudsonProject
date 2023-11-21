

selFile = 'U:\projects\2013_UnivMD_Maryland_71485\Analysis\NARW Det Eval\2023_Updated_Paths\MD01_Truth_logs\MD01_truth_selections.txt';
baseDir = 'U:\projects\2013_UnivMD_Maryland_71485\Sounds\71485_MD01_002K_10CH_AIFF';
thisTable = readtable(selFile,'Delimiter','\t','VariableNamingRule','preserve');

% For Windows
% dayFolders = dir(baseDir);
dataDays = regexp(thisTable.("Begin File"),'\d{8}','match');
dataDays = cellstr(char([dataDays{:}]));
dayFolders = strcat(repmat('002K_10CH_',length(dataDays),1),dataDays);
newPaths = fullfile(cellstr(repmat(baseDir,size(thisTable,1),1)),dayFolders,thisTable.("Begin File"));
saveStr = '_UpdatedPaths.txt';

% For Macs, also run this:
% newPaths = strrep(newPaths,'\','/'); 
% newPaths = strrep(newPaths,'U:','');
% newPaths = strcat(repmat('/Volumes/ag-clo-repnas5.ad.cornell.edu',length(newPaths),1),newPaths);
% saveStr = '_UpdatedPaths_Mac.txt';

thisTable.("Begin Path") = newPaths;
saveName = strrep(selFile,'.txt',saveStr);
writetable(thisTable,saveName,'Delimiter','\t');