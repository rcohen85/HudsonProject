% Test for correlation between fish counts and soundscape metrics

% Import fish count data
fishData20 = readtable('C:\Users\rec297\CCB\HudsonProject\BlackCreekHerringCountsPerHour.xlsx','Sheet',1);
fishData21 = readtable('C:\Users\rec297\CCB\HudsonProject\BlackCreekHerringCountsPerHour.xlsx','Sheet',2);

% calculate average fish/hr/day
countsPerHour20 = table2array(fishData20(:,11));
weirSets20 = table2array(fishData20(:,'daySet'))+days(table2array(fishData20(:,'timeSet')));
weirPulls20 = table2array(fishData20(:,'dayPull'))+days(table2array(fishData20(:,'timePull')));

% quick plot to check that data is continuous, no gaps (start of each set is end of previous set)
figure(1)
plot(weirSets20,ones(117,1),'o')
hold on
plot(weirPulls20,ones(117,1),'*')
hold off

weirDurs20 = hours(split(between(weirSets20,weirPulls20),'Time'));
weirDurs20(:,2) = floor(weirDurs20(:,1));
weirDurs20(:,3) = rem(weirDurs20(:,1),1);
allHours20 = (dateshift(weirSets20(1),'start','hour'):hours(1):dateshift(weirPulls20(end),'end','hour'))';
hourlyCounts20 = zeros(length(allHours20),1);
for i = 1:length(weirSets20)

    % find which hour bins overlap with this set
    thisChunk = find(allHours20>=dateshift(weirSets20(i),'start','hour') & allHours20<dateshift(weirPulls20(i),'end','hour'));

    if length(thisChunk)>=2
        % scale counts in first hour and last hours proportional to effort
        firstEff = hours(dateshift(weirSets20(i),'end','hour') - weirSets20(i));
        scaledFirstHrCount = countsPerHour20(i)*firstEff;
        lastEff = hours(weirPulls20(i) - dateshift(weirPulls20(i),'start','hour'));
        scaledLastHrCount = countsPerHour20(i)*lastEff;
        hourlyCounts20([thisChunk(1),thisChunk(end)],1) = [hourlyCounts20(thisChunk(1)) + scaledFirstHrCount,...
                                                            hourlyCounts20(thisChunk(end)) + scaledLastHrCount];

        if length(thisChunk)>=3
            % plug unscaled countsPerHour into all complete hours
            hourlyCounts20(thisChunk(2):thisChunk(end-1),1) = countsPerHour20(i);
        end

    elseif length(thisChunk)==1
        hourEff = hours(weirPulls20(i) - weirSets20(i));
        scaledCounts = countsPerHour20(i)*hourEff;
        hourlyCounts20(thisChunk) = hourlyCounts20(thisChunk) + scaledCounts;
    end
end

fishDates = [datetime(table2array(fishData20(:,1)),'InputFormat','d-MMM')];
fishCounts = (table2array(fishData20(:,15))./table2array(fishData20(:,14)))*24; % calculate hourly counts, then estimate daily counts

% Import broadband (BB) noise level data
BBdata = readtable('W:\projects\2022_NOAA-NERRS_HudsonNY_144488\SoundscapeMetrics\BlackCreek\BC05\BC05_BB_1h.csv');
BBdates = datetime(table2array(BBdata(:,1)),'InputFormat','yyyy-MM-dd''T''HH'); % dates get read in messed up from this one
BBlevels = str2double(extractAfter(table2array(BBdata(:,3)),'00.000Z,'));

% Import octave (O) levels
OLdata = readtable('W:\projects\2022_NOAA-NERRS_HudsonNY_144488\SoundscapeMetrics\BlackCreek\BC05\BC05_OL_1h.csv');
OLdates = datetime(table2array(OLdata(:,1)),'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z');
OLevels = table2array(OLdata(:,2:end));
OCenters = OLdata.Properties.VariableNames(2:end);
OCenters = str2double(strrep(cellfun(@(x) extractAfter(x,3),OCenters,'UniformOutput',0),'_','.'));

% Import third octave (TO) levels
TOLdata = readtable('W:\projects\2022_NOAA-NERRS_HudsonNY_144488\SoundscapeMetrics\BlackCreek\BC05\BC05_TOL_1h.csv');
TOLdates = datetime(table2array(TOLdata(:,1)),'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z');
TOLevels = table2array(TOLdata(:,2:end));
TOCenters = TOLdata.Properties.VariableNames(2:end);
TOCenters = str2double(strrep(cellfun(@(x) extractAfter(x,4),TOCenters,'UniformOutput',0),'_','.'));

% Compute daily average BB, O, and TO levels
dayBins = dateshift(BBdates(1),'start','day'):1:dateshift(BBdates(end),'end','day');
[N, edges,BBbin] = histcounts(BBdates,dayBins);
dailyBBlevels = grpstats(table(BBlevels,BBbin),'BBbin','mean');

[N, edges,OLbin] = histcounts(OLdates,dayBins);
dailyOlevels = grpstats(table(OLevels,OLbin),'OLbin','mean');

[N, edges,TOLbin] = histcounts(TOLdates,dayBins);
dailyTOlevels = grpstats(table(TOLevels,TOLbin),'TOLbin','mean');

% Bring in environmental data?

% Align dates and create matrix of counts and soundscape metrics
sharedDates = intersect(fishDates,BBdates);
fishInd = ismember(fishDates,sharedDates);
BBind = ismember(BBdates,sharedDates);
Oind = ismember(OLdates,sharedDates);
TOind = ismember(TOLdates,sharedDates);

dataMat = [fishCounts(fishInd),dailyBBlevels]

% Compute correlation matrix and plot
corrplot(dataMat)

