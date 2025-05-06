% Use annotations in selection tables to plot temporal occurrence of
% particular sounds of interest

%% Deployment-Specific Settings
selDir = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AtlStur\ManualSelectionTables'; % directory containing selection tables
% selDir = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AtlStur\ThunderClassifier_Output\CleanedForTimeseries';
selTabStr = '*ClackThunder.txt'; % wildcard to identify correct selection tables
% selTabStr = '*cleanedForTimeseries.txt'; % wildcard to identify correct selection tables
outDir = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AtlStur\Figures'; % directory to save figures
plotSite = 'Atlantic Sturgeon Spawning Site'; % site name for plot titles
saveSite = 'AtlSturSpwn'; % site name with no spaces for saved plot file names
labCol = 'Tags'; % column containing sound labels
label = 'Clack Thunder'; % sound label of interest in selection tables
site = [41.81440, -73.94540]; % lat, lon of recording location
% effort = datetime({'2021/11/19 11:28:17','2021/12/09 06:25:45'; % Shortnose sturgeon effort
%                     '2022/03/09 07:09:07','2022/03/29 01:51:05';
%                     '2022/10/12 13:17:04','2022/10/27 05:39:38'},'InputFormat','yyyy/MM/dd HH:mm:ss');
effort = datetime({'2021/04/24 09:05:03', '2021/04/25 21:04:05';
                    '2021/06/09 10:14:12',	'2021/06/22 13:27:53';
                    '2021/06/25 08:40:10',	'2021/06/30 20:37:52';
                    '2021/07/01 20:37:38','2021/07/12 23:35:44'},'InputFormat','yyyy/MM/dd HH:mm:ss'); % Atlantic sturgeon effort
tempRes = 'day'; % can be 'hour', 'day', or 'week' % temporal resolution for timeseries plot
UTCOffset = -4; % timezone for diel plot ***NOTE: IF EFFORT AND DETECTION DATA ARE ALREADY IN LOCAL TIME, SET UTCOFFSET TO 0 FOR DIEL PLOT

%% Compile calls of interest from selection tables
% load selections
selTabs = dir(fullfile(selDir,selTabStr));
detTimes = [];
for i=1:length(selTabs)
    % find annotations matching desired label
    tb = readtable(fullfile(selDir,selTabs(i).name),'Delimiter','\t','VariableNamingRule','preserve');
    specs = find(strcmp(table2array(tb(:,'View')),'Spectrogram 1'));
    tb = tb(specs,:);
    detInds = find(strcmp(table2array(tb(:,labCol)),label));

    % calculate absolute timestamps using filename timestamps and file offsets
    fileNames = table2array(tb(detInds,'Begin File'));
    fileStarts = cellfun(@(x) regexp(x,'\d{8}[_]\d{6}','match'),fileNames,'UniformOutput',0);
    fileStarts = (datetime([fileStarts{:}],'InputFormat','yyyyMMdd_HHmmss'))';
    detSt = fileStarts + seconds(table2array(tb(detInds,'File Offset (s)')));
    detEnd = detSt + seconds(table2array(tb(detInds,'End Time (s)'))-table2array(tb(detInds,'Begin Time (s)')));
    detTimes = [detTimes;detSt,detEnd];
    tb = [];
    specs = [];
    detInds = [];
    fileNames = [];
    fileStarts = [];
end

detTimes = sortrows(detTimes);

% Bin calls to 1-min intervals (to get around issue of detector not always counting individual calls, sometimes double-counting calls)
minBins = effort(1,1):minutes(1):effort(end,2);
[n,minBins,bin] = histcounts(detTimes(:,1),minBins);
minsWdets = minBins(find(n>0));

%% Binned timeseries plot
% bin detections to desire temporal resolution
binSt = dateshift(effort(1,1),'start',tempRes);
binEnd = dateshift(effort(end,2),'end',tempRes);
if strcmp(tempRes,'hour')
    int = hours(1);
elseif strcmp(tempRes,'day')
    int = days(1);
elseif strcmp(tempRes,'week')
    int = days(7);
end
timeBins = binSt:int:binEnd;
% [N,timeBins,bin] = histcounts(detTimes(:,1),timeBins);
% 
% % get rid of bins in off-effort periods
% if size(effort,1)==1
%     notOnEff = find(timeBins<dateshift(effort(1,1),'end',tempRes) | timeBins>dateshift(effort(1,2),'start',tempRes));
% else
%     for i = 1:(size(effort,1)-1)
%         gapBins = find(timeBins>=dateshift(effort(i,2),'start',tempRes) & timeBins<dateshift(effort(i+1,1),'end',tempRes));
%         timeBins(gapBins) = [];
%         N(gapBins) = [];
%     end
% end

% Bin call-positive-minutes to desired granularity
[minsPerInt,timeBins,bin] = histcounts(minsWdets,timeBins);

% Remove time bins outside of acoustic effort periods
badInds=[];
for i=1:length(timeBins)-1
    for j=1:size(effort,1)-1
        if timeBins(i)>effort(j,2) & timeBins(i)<effort(j+1,1)
            badInds = [badInds;i];
        end
    end
end
minsPerInt(badInds) = nan;

% Scale bins with partial effort


% plot
maxCount = max(minsPerInt);
if strcmp(tempRes,'day')
    period = 'Daily';
elseif strcmp(tempRes,'hour')
    period = 'Hourly';
elseif strcmp(tempRes,'week')
    period = 'Weekly';
end

figure(99), clf
    hold on
%     for i = 1:size(effort,1)
%         thisdep = find(timeBins>=dateshift(effort(i,1),'end',tempRes) & timeBins<dateshift(effort(i,2),'start',tempRes));
% %         if strcmp(tempRes,'hour')
% %         scatter(timeBins(thisdep),N(thisdep),'MarkerFaceColor','b','MarkerFaceAlpha',0.3,'MarkerEdgeColor','b','MarkerEdgeAlpha',0.3)
% %         else
%             plot(timeBins(thisdep),minsPerInt(thisdep),'-','LineWidth',2,'Color',[0 102 204]./255);
% %             plot(timeBins(thisdep),N(thisdep),'.','MarkerSize',2,'Color',[0 102 204]./255);
% %         end
%         
%     end
    plot(timeBins(1:end-1),minsPerInt,'-','LineWidth',2,'Color',[0 102 204]./255)
if size(effort,1)>1
    for i = 1:(size(effort,1)-1)
        noEff = patch([effort(i,2),effort(i,2),effort(i+1,1),effort(i+1,1)],[maxCount*1.1,0,0,maxCount*1.1],[128,128,128]./255);
        noEff.FaceAlpha = 0.1;
        noEff.EdgeColor = 'none';
    end
end
set(gca,'FontSize',14)
hold off
ylim([0, maxCount*1.1]);
xlim([effort(2,1),effort(end,2)]);
ylabel([period,' Count'])
title([label,'s at ',plotSite])

%save plot
saveName = [period,'_',label,'_Timeseries_at_',saveSite];
saveas(figure(99),fullfile(outDir,saveName),'png');
saveName = [saveName,'.pdf'];
exportgraphics(figure(99),fullfile(outDir,saveName),'ContentType','vector');

%% Diel plot
    [sunrise,sunset] = sunTimes(site(1,1),site(1,2),(effort(1,1)-days(1)):1:(effort(end,2)+days(1)),0);
    night = [sunset(1:end-1)',sunrise(2:end)'];
    effortGaps = [datetime,datetime];
    if size(effort,1)>1
        effortGaps(1,:) = [datetime('2020/01/01','InputFormat','yyyy/MM/dd'),effort(1,1)];
        for i = 1:(length(effort)-1)
            effortGaps(i+1,:) = [effort(i,2),effort(i+1,1)];
        end
        effortGaps(end+1,:) = [effort(i+1,2),datetime('2022/01/01','InputFormat','yyyy/MM/dd')];
    else
        effortGaps(1,:) = [datetime('2020/01/01','InputFormat','yyyy/MM/dd'),effort(1,1)];
        effortGaps(2,:) = [effort(1,2),datetime('2022/01/01','InputFormat','yyyy/MM/dd')];
    end
    

    figure(999), clf
    % add shading during nighttime hours
    [nightH,~,~] = visPresence(night, 'Color', 'black', ...
        'LineStyle', 'none', 'Transparency', .1, 'UTCOffset',UTCOffset,...
        'Resolution_m', 1, 'DateRange', [effort(2,1),effort(end,2)], 'DateTickInterval',5);
    

    % add shading during off-effort periods
    [effH,~,~] = visPresence(effortGaps, 'Color', 'purple', ...
        'LineStyle', 'none', 'Transparency', .1, 'UTCOffset',UTCOffset,...
        'Resolution_m', 1, 'DateRange', [effort(2,1),effort(end,2)], 'DateTickInterval',5);

    % add species presence data with solid bars
    [BarH, ~, ~] = visPresence(detTimes, 'Color','blue',...
        'UTCOffset',0,'Resolution_m',1, 'DateRange',[effort(2,1),effort(end,2)],...
        'DateTickInterval',3,'Title',[label,'s at ',plotSite]);
    set(gca,'FontSize',14)
%     set(gca,'YDir','reverse') % control whether earliest date is shown at the top or the bottom of the y-axis

    %save plot
    saveName = ['Diel_',label,'s_at_',saveSite];
    saveas(figure(999),fullfile(outDir,saveName),'png');
    saveName = [saveName,'.pdf'];
    exportgraphics(figure(999),fullfile(outDir,saveName),'ContentType','vector');


%% Diel histograms

[sunrise,sunset] = sunTimes(site(1,1),site(1,2),effort(1,1):1:(effort(end,2)+1),0); %sunrise/sunset in UTC
sunrise = sunrise + hours(UTCOffset);
sunset = sunset + hours(UTCOffset);
meanSR = mean(timeofday(sunrise));
meanSS = mean(timeofday(sunset));

% Unnormalized time of day
% % get time of day from full detection timestamps
% % ToD = datevec(detTimes(:,1));
% % ToD(:,1:3) = 0;
% % ToD = datetime(ToD);
% ToD2 = timeofday(detTimes(:,1));
% get time of day from call-positive-minute timestamps
ToD3 = timeofday(minsWdets);

% plot histogram of hourly bins with nighttime shaded
figure(22), clf
h = histogram(ToD3,24)
hold on
patch([0,meanSR,meanSR,0,0],[0,0,max(h.BinCounts)*1.1,max(h.BinCounts)*1.1,0],'k','FaceAlpha',0.2,'EdgeColor','none')
patch([meanSS,hours(24),hours(24),meanSS,meanSS],[0,0,max(h.BinCounts)*1.1,max(h.BinCounts)*1.1,0],'k','FaceAlpha',0.2,'EdgeColor','none')
title(['Hourly ',label,' Counts at ',plotSite])
ylabel('Counts');
xlim([hours(0),hours(24)]);
ylim([0,max(h.BinCounts)*1.1]);
xticks([hours(0),hours(6),hours(12),...
    hours(18),hours(24)])
xticklabels({'00:00','06:00','12:00','18:00','24:00'})
xlabel('Time of Day');
set(gca,'FontSize',14)
hold off

saveName = [label,'_DielHist_at_',saveSite];
saveas(figure(22),fullfile(outDir,saveName),'png');
saveName = [saveName,'.pdf'];
exportgraphics(figure(22),fullfile(outDir,saveName),'ContentType','vector');

% Normalized time of day, to account for day length changes
% normTimes = normTimeofD(detTimes(:,1),sunrise',sunset');
normTimes = normTimeofD(minsWdets,sunrise',sunset');

% plot histogram of hourly bins
figure(222),clf
h = histogram(normTimes,24);
hold on
% title(['Normalized ',label,' Times at ',plotSite])
patch([-0.02,1,1,-0.02,-0.02],[0,0,max(h.BinCounts)*1.1,max(h.BinCounts)*1.1,0],'k','FaceAlpha',0.2,'EdgeColor','none')
ylabel('Count of Calls');
title(['Hourly ',label,' Counts at ',plotSite])
ylim([0,max(h.BinCounts)*1.1]);
xlim([min(h.BinEdges),max(h.BinEdges)]);
xticks([-1,0,1]);
xticklabels({'Sunrise','Sunset','Sunrise'})
set(gca,'FontSize',14)
hold off

saveName = [label,'_NormDielHist_at_',saveSite];
saveas(figure(222),fullfile(outDir,saveName),'png');
saveName = [saveName,'.pdf'];
exportgraphics(figure(222),fullfile(outDir,saveName),'ContentType','vector');

%% Explore duration of calls - don't use when using CNN output detections w
% standardized durations, or when multiple calls may occur in a single detection

% % calculate durations
% detDurs = seconds(detTimes(:,2) - detTimes(:,1));
% 
% % histogram of call durations
% figure (33),clf
% histogram(detDurs,20)
% ylabel('Counts');
% xlabel('Duration (s)');
% title([label,' Durations at ',plotSite]);
% 
% saveName = [label,'_DurationHist_at_',saveSite];
% saveas(figure(33),fullfile(outDir,saveName),'png');
% 
% % plot durations as a function of time of day
% ToD = datevec(detTimes(:,1));
% ToD(:,1:3) = 0;
% ToD = datetime(ToD);
% 
% figure(333),clf
% scatter(ToD,detDurs,40,'MarkerFaceColor','b','MarkerFaceAlpha',0.3,'MarkerEdgeColor','b','MarkerEdgeAlpha',0.3)
% ylabel('Duration (s)');
% xlim([datetime(0,0,0,0,0,0),datetime(0,0,0,23,59,59)]);
% xticks([datetime(0,0,0,0,0,0),datetime(0,0,0,6,0,0),datetime(0,0,0,12,0,0),...
%     datetime(0,0,0,18,0,0),datetime(0,0,0,23,59,59)])
% xticklabels({'00:00','06:00','12:00','18:00','24:00'})
% xlabel('Time of Day');
% title([label,' Durations at ',plotSite]);
% saveName = [label,'_Duration_by_ToD_at_',saveSite];
% saveas(figure(333),fullfile(outDir,saveName),'png');


