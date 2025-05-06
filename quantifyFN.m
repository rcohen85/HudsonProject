% Quantify FN/hr to evaluate detector performance over time

selDir = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AtlStur\ThunderClassifier_Output\ErrorEval';
col = 'Score'; % column containing flag for FNs
label = 'Manual'; % label identifying FNs
effort = datetime({'2021/06/09 10:14:12',	'2021/06/22 13:27:53'; % Period(s) of time which were reviewed for FN 
                    '2021/06/25 08:40:10',	'2021/06/30 20:37:52';
                    '2021/07/01 20:37:38','2021/07/12 23:35:44'},'InputFormat','yyyy/MM/dd HH:mm:ss');

selTabs = dir(fullfile(selDir,'*FN.txt'));

FNtimes = [];

% Aggregate timestamps of FN detections across selection tables
for i=1:length(selTabs)
    tb = readtable(fullfile(selDir,selTabs(i).name),'Delimiter','\t','VariableNamingRule','preserve');
    specs = find(strcmp(table2array(tb(:,'View')),'Spectrogram 1'));
    tb = tb(specs,:);

    FNinds = find(strcmp(table2array(tb(:,col)),label));

    fileNames = table2array(tb(FNinds,'Begin File'));
    fileStarts = cellfun(@(x) regexp(x,'\d{8}[_]\d{6}','match'),fileNames,'UniformOutput',0);
    fileStarts = (datetime([fileStarts{:}],'InputFormat','yyyyMMdd_HHmmss'))';
    detSt = fileStarts + seconds(table2array(tb(FNinds,'File Offset (s)')));

    FNtimes = [FNtimes;detSt];

end

% Bin FN to hourly resolution
hourBins = effort(1,1):hours(1):effort(end,2);
[hourlyFN,hourBins,bin] = histcounts(FNtimes(:,1),hourBins);

% Remove time bins outside of acoustic effort periods
badInds=[];
for i=1:length(hourBins)-1
    for j=1:size(effort,1)-1
        if hourBins(i)>effort(j,2) & hourBins(i)<effort(j+1,1);
            badInds = [badInds;i];
        end
    end
end
hourlyFN(badInds) = nan;
maxCount = max(hourlyFN);

figure(111)
plot(hourBins(1:end-1),hourlyFN,'o')
hold on
if size(effort,1)>1
    for i = 1:(size(effort,1)-1)
        noEff = patch([effort(i,2),effort(i,2),effort(i+1,1),effort(i+1,1)],[maxCount*1.1,0,0,maxCount*1.1],[128,128,128]./255);
        noEff.FaceAlpha = 0.1;
        noEff.EdgeColor = 'none';
    end
end
hold off
ylim([0, maxCount*1.1]);
xlim([effort(1,1),effort(end,2)]);
ylabel('Hourly FN Count');
set(gca,'FontSize',14)
