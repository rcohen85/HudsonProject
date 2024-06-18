
selDir = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AtlStur\ThunderClassifier_Output';
fileNameWildcard = '*cleanedForTimeseries';
label = 'Thunder'; % label of interest
labCol = 'Tags'; % column containing sound labels
telemData = readtable('W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AtlStur\Telemetry_Counts_EsopusIsland_2021.xlsx','VariableNamingRule','preserve');
timeInt = 0.5; % temporal granularity of telemetry data (hours)
newInt = 24; % if desired to re-bin data to lower temporal granularity, enter desired resolution (hours); otherwise, leave empty

% wrangle telemetry data
telemTimes = table2array(telemData(:,'Date.and.Time..UTC.'));
telemCounts = table2array(telemData(:,'Number of unique detections'));
if isempty(newInt)
fullTimeBins = telemTimes(1):hours(timeInt):(telemTimes(end)+hours(timeInt));
[lia, locb] = ismember(telemTimes,fullTimeBins);
fullTelemCounts = zeros(length(fullTimeBins)-1,1);
fullTelemCounts(locb) = telemCounts(lia);
else
fullTimeBins = telemTimes(1):hours(newInt):(telemTimes(end)+hours(newInt));
[n,fullTimeBins,bin] = histcounts(telemTimes,fullTimeBins);
[G,ID] = findgroups(bin);
newTelemCounts = splitapply(@sum,telemCounts,G);
newTelemTimes = fullTimeBins(unique(bin));
[lia, locb] = ismember(newTelemTimes,fullTimeBins);
fullTelemCounts = zeros(length(fullTimeBins)-1,1);
fullTelemCounts(locb) = newTelemCounts(lia);
end

% Load call occurrence data
selList = dir(fullfile(selDir,[fileNameWildcard,'.txt']));
detStarts = [];
detEnds = [];
for i=1:length(selList)
    tb = readtable(fullfile(selDir,selList(i).name),'Delimiter','\t','VariableNamingRule','preserve');
    specs = find(strcmp(table2array(tb(:,'View')),'Spectrogram 1'));
    tb = tb(specs,:);
    detInds = find(strcmp(table2array(tb(:,labCol)),label));

    fileNames = table2array(tb(detInds,'Begin Path'));
    fileStarts = cellfun(@(x) regexp(x,'\d{8}[_]\d{6}','match'),fileNames,'UniformOutput',0);
    fileStarts = (datetime([fileStarts{:}],'InputFormat','yyyyMMdd_HHmmss'))';
    detSt = fileStarts + seconds(table2array(tb(detInds,'File Offset (s)')));
    detEnd = detSt + seconds(table2array(tb(detInds,'End Time (s)'))-table2array(tb(detInds,'Begin Time (s)')));
    detStarts = [detStarts;detSt];
    detEnds = [detEnds;detEnd];

end

% % Cut detections up to max 1 min
% longDets = find(detEnds>=dateshift(detStarts,'end','minute'));
% if ~isempty(longDets)
%     q = dateshift(detStarts,'end','minute');
%     if longDets(1)==1
%         starts = [detStarts(1);q(1)];
%         ends = [q(longDets(1));detEnds(1)];
%     else
%         starts = detStarts(1:longDets(1));
%         ends = detEnds(1:longDets(1)-1);
%         starts = [starts;q(longDets(1))];
%         ends = [ends;q(longDets(1));detEnds(longDets(1))];
% 
%     end
%     if length(longDets)>1
%         for i=2:length(longDets)
%             starts = [starts;q(longDets(i))];
%             ends = [ends;q(longDets(i));detEnds(longDets(i))];
%         end
%     else
%     end
% 
% else
%     starts = detStarts;
%     ends = detEnds;
% end

% Bin calls to 1-min intervals
minBins = telemTimes(1):minutes(1):(telemTimes(end)+minutes(1));
[n,minBins,bin] = histcounts(detStarts,minBins);
minsWdets = minBins(find(n>0));

% Bin call-positive-minutes to match granularity of telemetry data
[minsPerInt,fullTimeBins,bin] = histcounts(minsWdets,fullTimeBins);

% % Bin calls to match granularity of telemetry data
% [detsPerInt,fullTimeBins,bin] = histcounts(detStarts,fullTimeBins);

% Calculate correlation
rho = corr(fullTelemCounts,movmean(minsPerInt,10)');
% fprintf('Rho = %0.2f\n',rho);

% Plot telemetry counts and call-positive-minutes per interval
figure(99)
plot(fullTimeBins(1:end-1),fullTelemCounts)
hold on
plot(fullTimeBins(1:end-1),minsPerInt)
% plot(fullTimeBins(1:end-1),detsPerInt)
ylabel('Counts')
if isempty(newInt)
    title([num2str(timeInt),'-hour Bin'])
else
    title([num2str(newInt),'-hour Bin'])
end
legend('# Unique Telemetry Tags','# Call-positive-minutes')
text(fullTimeBins(round(length(fullTimeBins)*0.85)),max([fullTelemCounts;minsPerInt'])*0.8,sprintf('Rho = %0.2f\n',rho),'FontSize',14)
set(gca,'FontSize',14)
hold off


