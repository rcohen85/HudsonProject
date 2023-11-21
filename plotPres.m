% Use annotations in selection tables to plot temporal occurrence of
% particular sounds of interest

%% Deployment-Specific Settings
selDir = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AtlStur\SelectionTables'; % directory containing selection tables
selTabStr = '*ManualSelections.txt'; % wildcard to identify correct selection tables
outDir = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AtlStur\Figures'; % directory to save figures
plotSite = 'Atlantic Sturgeon Spawning Site'; % site name for plot titles
saveSite = 'AtlSturSpwn'; % site name with no spaces for saved plot file names
labCol = 'SoundType'; % column containing sound labels
label = 'Rumble'; % sound label of interest in selection tables
site = [41.81440, -73.94540]; % lat, lon of recording location
% effort = datetime({'2021/11/19 11:28:17','2021/12/09 06:25:45'; % Shortnose sturgeon effort
%                     '2022/03/09 07:09:07','2022/03/29 01:51:05';
%                     '2022/10/12 13:17:04','2022/10/27 05:39:38'},'InputFormat','yyyy/MM/dd HH:mm:ss');
effort = datetime({'2021/06/09 10:14:12',	'2021/06/22 13:27:53'; % Atlantic sturgeon effort
                    '2021/06/25 08:40:10',	'2021/06/30 20:37:52';
                    '2021/07/01 20:37:38',	'2021/07/12 20:35:44';
                    '2022/05/27 06:49:50',	'2022/05/28 06:48:36';
                    '2022/06/16 06:44:28',	'2022/06/16 07:36:31';
                    '2022/08/03 11:09:30',	'2022/08/10 10:51:27'},'InputFormat','yyyy/MM/dd HH:mm:ss');
tempRes = 'hour'; % can be 'hour', 'day', or 'week' % temporal resolution for timeseries plot
UTCOffset = -4; % timezone for diel plot

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

%% Binned timeseries plot
% bin detections to desire temporal resolution (excluding partial effort bins)
binSt = dateshift(effort(1,1),'end',tempRes);
binEnd = dateshift(effort(end,2),'start',tempRes);
if strcmp(tempRes,'hour')
    timeBins = binSt:hours(1):binEnd;
elseif strcmp(tempRes,'day')
    timeBins = binSt:days(1):binEnd;
elseif strcmp(tempRes,'week')
    timeBins = binSt:days(7):binEnd;
end
[N,timeBins,bin] = histcounts(detTimes(:,1),timeBins);

% get rid of bins in off-effort periods
for i = 1:(length(effort)-1)
gapBins = find(timeBins>=dateshift(effort(i,2),'start',tempRes) & timeBins<dateshift(effort(i+1,1),'end',tempRes));
timeBins(gapBins) = [];
N(gapBins) = [];
end

% plot
maxCount = max(N);
if strcmp(tempRes,'day')
    period = 'Daily';
elseif strcmp(tempRes,'hour')
    period = 'Hourly';
elseif strcmp(tempRes,'week')
    period = 'Weekly';
end

figure(99), clf
hold on
for i = 1:length(effort)
    thisdep = find(timeBins>=dateshift(effort(i,1),'end',tempRes) & timeBins<dateshift(effort(i,2),'start',tempRes));
    if any(strcmp(tempRes,{'day','week'}))
    plot(timeBins(thisdep),N(thisdep),'-','LineWidth',2,'Color',[0 102 204]./255);
    elseif strcmp(tempRes,'hour')
    scatter(timeBins(thisdep),N(thisdep),'MarkerFaceColor','b','MarkerFaceAlpha',0.3,'MarkerEdgeColor','b','MarkerEdgeAlpha',0.3)
    end
end
for i = 1:(length(effort)-1)
    noEff = patch([effort(i,2),effort(i,2),effort(i+1,1),effort(i+1,1)],[maxCount*1.1,0,0,maxCount*1.1],[128,128,128]./255);
    noEff.FaceAlpha = 0.2;
    noEff.EdgeColor = 'none';
end
hold off
ylim([0, maxCount*1.1]);
xlim([timeBins(1),timeBins(end-1)]);
ylabel([period,' Count'])
title([label,'s at ',plotSite])

%save plot
saveName = [period,'_',label,'_Timeseries_at_',saveSite,'_3'];
saveas(figure(99),fullfile(outDir,saveName),'png');

%% Diel plot
    [sunrise,sunset] = sunTimes(site(1,1),site(1,2),effort(1,1):1:effort(end,2),0);
    night = [sunset(1:end-1)',sunrise(2:end)'];
    effortGaps = [datetime,datetime];
    for i = 1:(length(effort)-1)
        effortGaps(i,:) = [effort(i,2),effort(i+1,1)];
    end

    figure(999), clf
    % add shading during nighttime hours
    [nightH,~,~] = visPresence(night, 'Color', 'black', ...
        'LineStyle', 'none', 'Transparency', .15, 'UTCOffset',UTCOffset,...
        'Resolution_m', 1, 'DateRange', [effort(1,1),effort(end,2)], 'DateTickInterval',90);

    % add shading during off-effort periods
    [effH,~,~] = visPresence(effortGaps, 'Color', 'purple', ...
        'LineStyle', 'none', 'Transparency', .15, 'UTCOffset',UTCOffset,...
        'Resolution_m', 1, 'DateRange', [effort(1,1),effort(end,2)], 'DateTickInterval',90);

    % add species presence data with solid bars
    [BarH, ~, ~] = visPresence(detTimes, 'Color','blue',...
        'UTCOffset',UTCOffset,'Resolution_m',5, 'DateRange',[effort(1,1),effort(end,2)],...
        'DateTickInterval',60,'Title',[label,'s at ',plotSite]);

    %save plot
    saveName = ['Diel_',label,'s_at_',saveSite];
    saveas(figure(999),fullfile(outDir,saveName),'png');


%% Diel histograms

% Unnormalized time of day
% get time of day from full detection timestamps
ToD = datevec(detTimes(:,1));
ToD(:,1:3) = 0;
ToD = datetime(ToD);

% plot histogram of hourly bins
figure(22)
histogram(ToD,24)
title(['Diel ',label,' Counts at ',plotSite])
ylabel('Counts');
xlim([datetime(0,0,0,0,0,0),datetime(0,0,0,23,59,59)]);
xticks([datetime(0,0,0,0,0,0),datetime(0,0,0,6,0,0),datetime(0,0,0,12,0,0),...
    datetime(0,0,0,18,0,0),datetime(0,0,0,23,59,59)])
xticklabels({'00:00','06:00','12:00','18:00','24:00'})
xlabel('Time of Day');

saveName = [label,'_DielHist_at_',saveSite];
saveas(figure(22),fullfile(outDir,saveName),'png');


% Normalized time of day, to account for day length changes
[sunrise,sunset] = sunTimes(site(1,1),site(1,2),effort(1,1):1:(effort(end,2)+1),0); %sunrise/sunset in UTC
sunrise = sunrise + hours(UTCOffset);
sunset = sunset + hours(UTCOffset);
normTimes = normTimeofD(detTimes(:,1),sunrise',sunset');

% plot histogram of hourly bins
figure(222)
histogram(normTimes,24)
title(['Normalized ',label,' Times at ',plotSite])
ylabel('Counts');
xlim([-1,1]);
xticks([-1,0,1]);
xticklabels({'Sunrise','Sunset','Sunrise'})

saveName = [label,'_NormDielHist_at_',saveSite];
saveas(figure(222),fullfile(outDir,saveName),'png');


%% Explore duration of calls

% calculate durations
detDurs = seconds(detTimes(:,2) - detTimes(:,1));

% histogram of call durations
figure (33),clf
histogram(detDurs,20)
ylabel('Counts');
xlabel('Duration (s)');
title([label,' Durations at ',plotSite]);

saveName = [label,'_DurationHist_at_',saveSite];
saveas(figure(33),fullfile(outDir,saveName),'png');

% plot durations as a function of time of day
ToD = datevec(detTimes(:,1));
ToD(:,1:3) = 0;
ToD = datetime(ToD);

figure(333),clf
scatter(ToD,detDurs,40,'MarkerFaceColor','b','MarkerFaceAlpha',0.3,'MarkerEdgeColor','b','MarkerEdgeAlpha',0.3)
ylabel('Duration (s)');
xlim([datetime(0,0,0,0,0,0),datetime(0,0,0,23,59,59)]);
xticks([datetime(0,0,0,0,0,0),datetime(0,0,0,6,0,0),datetime(0,0,0,12,0,0),...
    datetime(0,0,0,18,0,0),datetime(0,0,0,23,59,59)])
xticklabels({'00:00','06:00','12:00','18:00','24:00'})
xlabel('Time of Day');
title([label,' Durations at ',plotSite]);
saveName = [label,'_Duration_by_ToD_at_',saveSite];
saveas(figure(333),fullfile(outDir,saveName),'png');


