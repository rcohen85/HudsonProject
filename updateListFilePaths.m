clearvars
opts = delimitedTextImportOptions('DataLines',1,'Delimiter','tab');
inFile = readtable('R:\projects\2020_UnivMD_Maryland_139635\Williamson_HonorsThesis\OtherData\2015_BOEMVA_76085\76085_BOEMVA03_002K_FLAC_listfile.txt',opts);
saveName = 'R:\projects\2020_UnivMD_Maryland_139635\Williamson_HonorsThesis\OtherData\2015_BOEMVA_76085\76085_BOEMVA03_002K_FLAC_listfile_Mac.txt';

% newPaths = strrep(inFile.Var1,'W:\projects\2020_UnivMD_Maryland_139635\01_139635_FLAC\139635_Dep01\005K\','W:\projects\2020_UnivMD_Maryland_139635\139635_Dep01\01_139635_FLAC\139635_Dep01_005K_FLAC\');
newPaths = strrep(inFile.Var1,'\','/');
newPaths = strrep(newPaths,'R:','');
newPaths = strcat(repmat('/Volumes/ag-clo-repnas5.ad.cornell.edu',length(newPaths),1),newPaths);
newTable = cell2table(newPaths);

writetable(newTable,saveName,'Delimiter','\t','WriteVariableNames',0);



%% Make list files for selected days

inDir = 'R:\projects\2020_UnivMD_Maryland_139635\Williamson_HonorsThesis\OtherData\2011_NEaq_MA-RI_65895\01_Sounds\65895_MCEC03_20130212\CH01';
saveName = 'R:\projects\2020_UnivMD_Maryland_139635\Williamson_HonorsThesis\OtherData\2011_NEaq_MA-RI_65895\NewListFiles\2011_NEaq_MA-RI_65895_CH01_002K_FLAC_listfile.txt';
% annotDates = {'20141103','20141113','20141123','20141213','20141213','20141223','20150102','20150112','20150122'... % MD 01
%     '20150201','20150211','20150221','20150302','20150312','20150322','20150331','20150409','20150419'};
% annotDates = {'20160225','20160306','20160326','20160405','20160415','20160505',... % MD04
%     '20160515','20160525','20160604','20160614','20160624','20160704','20160714'};
annotDates = [datetime(2013,02,16):days(1):datetime(2013,07,30)];
annotDates = cellstr(datestr(annotDates,'yyyymmdd'));

dayDirs = dir(inDir);
dirNames = cellstr(vertcat(dayDirs(3:end).name));
fileNames = [];

for i = 1:length(annotDates)

    dirInd = find(contains(dirNames,annotDates{i}));
    theseFiles = dir(fullfile(inDir,dirNames{dirInd},'**/*.flac'));
    fullPaths = fullfile(cellstr(vertcat(theseFiles.folder)),cellstr(vertcat(theseFiles.name)));
    fileNames = [fileNames;fullPaths];

end


listTable = cell2table(fileNames);
writetable(listTable,saveName,'Delimiter','\t','WriteVariableNames',0)



%% Make list files for all days
clearvars
inDir = 'R:\projects\2020_UnivMD_Maryland_139635\Williamson_HonorsThesis\OtherData\2011_NEaq_MA-RI_65895\01_Sounds\65895_MCEC06_20140923\CH01';
saveName = 'R:\projects\2020_UnivMD_Maryland_139635\Williamson_HonorsThesis\OtherData\2011_NEaq_MA-RI_65895\NewListFiles\2011_NEaq_MA-RI_65895_MCEC06_CH01_002K_FLAC_listfile.txt';

allFiles = dir(fullfile(inDir,'**/*.flac'));
fullPaths = fullfile(cellstr(vertcat(allFiles.folder)),cellstr(vertcat(allFiles.name)));

listTable = cell2table(fullPaths);
writetable(listTable,saveName,'Delimiter','\t','WriteVariableNames',0)