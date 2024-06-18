%% Convert list files from Windows path format to Mac path format

opts = delimitedTextImportOptions('DataLines',1,'Delimiter','tab');
inFile = readtable('W:\projects\2020_UnivMD_Maryland_139635\03_139635_Analysis\Dep01\Species_Analyses\ManualAnalysisListFiles\T3M_RH419_ManualAnalysis.txt',opts);
saveName = 'W:\projects\2020_UnivMD_Maryland_139635\03_139635_Analysis\Dep01\Species_Analyses\ManualAnalysisListFiles\T3M_RH419_ManualAnalysis_ClaireMac.txt';

newPaths = strrep(inFile.Var1,'\','/');
% newPaths = strrep(newPaths,'W:','');
% newPaths = strcat(repmat('/Volumes/ag-clo-repnas5.ad.cornell.edu',length(newPaths),1),newPaths);
slashInd = strfind(newPaths,'/');
slashInd = cell2mat(cellfun(@(x) x(end),slashInd,'UniformOutput',0));
newPaths = strcat(repmat('/Volumes/Extreme Pro/2022_UnivMD_Maryland_139635/Sounds/T3M_RH419/',length(newPaths),1),extractAfter(newPaths,slashInd));
newTable = cell2table(newPaths);

writetable(newTable,saveName,'Delimiter','\t','WriteVariableNames',0);



%% Make list files for selected days

inDir = 'W:\projects\2020_UnivMD_Maryland_139635\01_139635_FLAC\139635_Dep01\005K\T3M_RH419';
saveName = 'W:\projects\2020_UnivMD_Maryland_139635\03_139635_Analysis\Dep01\Species_Analyses\ManualAnalysisListFiles\T3M_RH419_ManualAnalysis.txt';
fileExt = '.flac';
% annotDates = {'20141103','20141113','20141123','20141213','20141213','20141223','20150102','20150112','20150122'... % MD 01
%     '20150201','20150211','20150221','20150302','20150312','20150322','20150331','20150409','20150419'};
% annotDates = {'20160225','20160306','20160326','20160405','20160415','20160505',... % MD04
%     '20160515','20160525','20160604','20160614','20160624','20160704','20160714'};
% annotDates = {'20170612','20170622','20170702','20170712','20170722','20170801',...
%     '20170811','20170821','20170831','20170910','20170920','20170930','20171001',...
%     '20171011','20171021','20171031','20171110','20171120','20171130','20171210',...
%     '20171220','20171230','20180109'};
dateInts = {'20230620','20230831';'20230901','20231009'};
pct = [0.1,0.25];
annotDates = [];
for i = 1:size(dateInts,1)
    tot = days(datetime(dateInts{i,2},'InputFormat','yyyyMMdd')-datetime(dateInts{i,1},'InputFormat','yyyyMMdd'));
    int = round(tot/round(tot*pct(i)));
    annotDates = [annotDates;(datetime(dateInts{i,1},'InputFormat','yyyyMMdd'):int:datetime(dateInts{i,2},'InputFormat','yyyyMMdd'))'];
end
annotDates = cellstr(datestr(annotDates,'yyyymmdd'));

dayDirs = dir(inDir);
dirNames = cellstr(vertcat(dayDirs(3:end).name));
fileNames = [];

for i = 1:length(annotDates)

    dirInd = find(contains(dirNames,annotDates{i}));
    theseFiles = dir(fullfile(inDir,dirNames{dirInd},['*' fileExt]));
    fullPaths = fullfile(cellstr(vertcat(theseFiles.folder)),cellstr(vertcat(theseFiles.name)));
    fileNames = [fileNames;fullPaths];

end


listTable = cell2table(fileNames);
writetable(listTable,saveName,'Delimiter','\t','WriteVariableNames',0)



%% Make list files for all days

inDir = 'U:\projects\2013_UnivMD_Maryland_71485\Sounds\71485_MD07_002K_12CH_AIFF';
saveName = 'U:\projects\2013_UnivMD_Maryland_71485\KooguNARWDetEval\DetectorListFiles\MD07_listfile.txt';

allFiles = dir(fullfile(inDir,'**/*.aif'));
fullPaths = fullfile(cellstr(vertcat(allFiles.folder)),cellstr(vertcat(allFiles.name)));

listTable = cell2table(fullPaths);
writetable(listTable,saveName,'Delimiter','\t','WriteVariableNames',0)