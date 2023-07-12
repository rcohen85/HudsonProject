% Remove spurious automatic detections during periods of boat/mechanical noise

% directory containing pulse detector output
detDir = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\ShortStur\PulseDetectorOutputmetadata\SS02';
% selection table containing boat/noise times
man = readtable('W:\projects\2022_NOAA-NERRS_HudsonNY_144488\ShortStur\SelectionTables\144488Hudson_024K_SS02.ManualSelections.txt','Delimiter','tab');

%%
cleanDir = fullfile(detDir,'CleanedDetections');
if ~isdir(cleanDir)
    mkdir(cleanDir)
end

% Remove duplicate entries from boat/noise annotations (due to multiple sound views in Raven)
manWav = find(strcmp(table2array(man(:,'View')),'Spectrogram 1'));
man = man(manWav,:);

noiseInds = find(strcmp(table2array(man(:,'SoundType')),'Boat')|...
    strcmp(table2array(man(:,'SoundType')),'Mechanical Noise')|...
    strcmp(table2array(man(:,'SoundType')),'Deployment Noise'));
noiseFiles = table2array(man(noiseInds,'BeginFile'));
noiseStarts = table2array(man(noiseInds,'FileOffset_s_'));
noiseEnds = noiseStarts + table2array(man(noiseInds,'DeltaTime_s_'));

% loop through detection files and create new ones without detections
% during boat/noise events
detFiles = dir(fullfile(detDir,'*.mat'));

for i = 1:length(detFiles) % loop through detection files
    thisFile = detFiles(i).name;
    noiseEvents = find(strcmp(noiseFiles,strrep(thisFile,'.mat','.flac'))); % find noise events in each file
    load(fullfile(detDir,thisFile));

    if ~isempty(noiseEvents)
        for j = 1:length(noiseEvents) % identify and remove detections during noise
            badInds = find(clickTimes(:,1)>=noiseStarts(noiseEvents(j)) & clickTimes(:,1)<=noiseEnds(noiseEvents(j)));
            keepInds = setdiff(1:size(clickTimes,1),badInds);
            clickTimes = clickTimes(keepInds,:);
            ppSignal = ppSignal(keepInds,:);
            specClickTf = specClickTf(keepInds,:);
            yFiltBuff = yFiltBuff(keepInds,:);
        end
    end

    % save new, cleaned detection file
    save(fullfile(cleanDir,thisFile),'clickTimes','ppSignal','f','hdr','specClickTf',...
        'yFiltBuff','p','-v7.3');

end