
selDir = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AtlStur\Sound_Exploration_Characterization\ThunderClassifier_Output\CleanedForTimeseries';
fileNameWildcard = '*cleanedForTimeseries';
label = 'Thunder'; % label of interest
labCol = 'Tags'; % column containing sound labels
telemData = readtable("W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AtlStur\Sound_Exploration_Characterization\TelemetryDetections_EsopusIsland_2021.xlsx",'VariableNamingRule','preserve');
timeInt = 24; % bin data to desired temporal granularity(hours)
normalize = 1; % normalize timeseries for plotting?
acousticEffort = datetime({'2021/06/10','2021/06/22';
                            '2021/06/26','2021/06/30';
                            '2021/07/02','2021/07/12'},'InputFormat','uuuu/MM/dd');
saveDir = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AtlStur\Sound_Exploration_Characterization\Figures\Revision';

% wrangle telemetry data
telemTimes = table2array(telemData(:,'Date and Time (UTC)'));
telemPinger = table2array(telemData(:,'Transmitter'));
% telemPinger = str2double(extractAfter(telemPinger,'A69-9001-'));
fullTimeBins = dateshift(telemTimes(1),'start','day')+1:hours(timeInt):dateshift(telemTimes(end),'start','day');
[n,fullTimeBins,bin] = histcounts(telemTimes,fullTimeBins);
telemTimes(bin==0) = [];
telemPinger(bin==0) = [];
bin(bin==0) = [];
[G,ID] = findgroups(bin);
uniqueCounts = @(x) size(unique(x),1);
newTelemCounts = splitapply(uniqueCounts,telemPinger,G);
newTelemTimes = fullTimeBins(unique(bin));
[lia, locb] = ismember(newTelemTimes,fullTimeBins);
fullTelemCounts = zeros(length(fullTimeBins)-1,1);
fullTelemCounts(locb) = newTelemCounts(lia);

% Load call occurrence data
selList = dir(fullfile(selDir,[fileNameWildcard,'.txt']));
detStarts = [];
detEnds = [];
for i=1:length(selList)
    tb = readtable(fullfile(selDir,selList(i).name),'Delimiter','\t','VariableNamingRule','preserve');
    specs = find(strcmp(table2array(tb(:,'View')),'Spectrogram 1'));
    tb = tb(specs,:);
    detInds = find(strcmp(table2array(tb(:,labCol)),label));

    fileNames = table2array(tb(detInds,'Begin File'));
    fileStarts = cellfun(@(x) regexp(x,'\d{8}[_]\d{6}','match'),fileNames,'UniformOutput',0);
    fileStarts = (datetime([fileStarts{:}],'InputFormat','yyyyMMdd_HHmmss'))';
    detSt = fileStarts + seconds(table2array(tb(detInds,'File Offset (s)')));
    detEnd = detSt + seconds(table2array(tb(detInds,'End Time (s)'))-table2array(tb(detInds,'Begin Time (s)')));
    detStarts = [detStarts;detSt];
    detEnds = [detEnds;detEnd];

end

% Bin calls to 1-min intervals
minBins = telemTimes(1):minutes(1):(telemTimes(end)+minutes(1));
[n,minBins,bin] = histcounts(detStarts,minBins);
minsWdets = minBins(find(n>0));

% Bin call-positive-minutes to match granularity of telemetry data
[minsPerInt,fullTimeBins,bin] = histcounts(minsWdets,fullTimeBins);

% Remove time bins outside of acoustic effort periods
badInds=[];
for i=1:length(fullTimeBins)-1
    for j=1:size(acousticEffort,1)-1
        if fullTimeBins(i)>acousticEffort(j,2) & fullTimeBins(i)<acousticEffort(j+1,1)
            badInds = [badInds;i];
        end
    end
end
minsPerInt(badInds) = nan;

% Calculate correlation
% rho = corrcoef(fullTelemCounts,movmean(minsPerInt,10)','Rows','complete');
[rho,p] = corrcoef(fullTelemCounts,minsPerInt,'Rows','complete');
rho = rho(1,2);
p = p(1,2);

% Normalize Telem counts & Thunder counts for better visualization
if normalize
    norm_TelemCounts = fullTelemCounts-min(fullTelemCounts);
    norm_TelemCounts = norm_TelemCounts./max(norm_TelemCounts);
    norm_minsPerInt = minsPerInt-min(minsPerInt);
    norm_minsPerInt = norm_minsPerInt./max(norm_minsPerInt);

    % Plot telemetry counts and call-positive-minutes per interval
    figure(999),clf
    hold on
    for i=1:size(acousticEffort,1)-1
        patch([acousticEffort(i,2),acousticEffort(i,2),acousticEffort(i+1,1), ...
            acousticEffort(i+1,1),acousticEffort(i,2)],[0,1,1,0,0],'k','FaceAlpha',0.075,'EdgeColor','none')
    end
    plot(fullTimeBins(1:end-1),norm_TelemCounts,'LineWidth',2) 
    plot(fullTimeBins(1:end-1),norm_minsPerInt,'LineWidth',2)
    ylabel('Normalized Counts')
    %title([num2str(timeInt),'-hour Bin'])
    legend('','','Fish Abundance Index','Thunder-Positive Minutes')
    text(fullTimeBins(round(length(fullTimeBins)*0.83)),0.65,sprintf('Rho = %0.2f\n',rho),'FontSize',14)
    set(gca,'FontSize',14)
    hold off
    xlim([fullTimeBins(1),fullTimeBins(end-1)])

    savename = fullfile(saveDir,[strrep(label,' ',''),'TelemCorr_Daily_Norm.png']);
    saveas(figure(999),savename)
    savename = fullfile(saveDir,[strrep(label,' ',''),'TelemCorr_Daily_Norm.pdf']);
    exportgraphics(figure(999),savename,'ContentType','vector');

else
    figure(999),clf
     hold on
    for i=1:size(acousticEffort,1)-1
        patch([acousticEffort(i,2),acousticEffort(i,2),acousticEffort(i+1,1), ...
            acousticEffort(i+1,1),acousticEffort(i,2)],[0,max([fullTelemCounts;minsPerInt']),max([fullTelemCounts;minsPerInt']),0,0], ...
            'k','FaceAlpha',0.075,'EdgeColor','none')
    end
    plot(fullTimeBins(1:end-1),fullTelemCounts,'LineWidth',2)
    plot(fullTimeBins(1:end-1),minsPerInt,'LineWidth',2)
    ylabel('Counts')
    title([num2str(timeInt),'-hour Bin'])
    legend('','','# Unique Telemetry Tags','# Call-positive-minutes')
    text(fullTimeBins(round(length(fullTimeBins)*0.85)),max([fullTelemCounts;minsPerInt'])*0.65, ...
        sprintf('Rho = %0.2f\n',rho),'FontSize',14)
    set(gca,'FontSize',14)
    hold off
    xlim([fullTimeBins(1),fullTimeBins(end-1)])
    ylim([0,max([fullTelemCounts;minsPerInt'])])

    savename = fullfile(saveDir,[strrep(label,' ',''),'TelemCorr_Daily.png']);
    saveas(figure(999),savename)
    savename = fullfile(saveDir,[strrep(label,' ',''),'TelemCorr_Daily.pdf']);
    exportgraphics(figure(999),savename,'ContentType','vector');
end

