% Create short clips with annotated signals centered in each
% Works best if annotations are shorter than/similar to the length of the
% clip (winSize); if winSize is set shorter than some annotation durations,
% those annotations will be truncated. Separate folders will be created for
% each unique label in the selection table(s)

% Path to audio files
audioDir = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_FLAC\AtlSturSpw\ST300'; 
% file extension of audio files
fileExt = '.flac'; 
% Path to selection table(s) with annotations to be clipped
selDir = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AtlStur\ThunderClassifier_Output\CleanedForTimeseries\Noise'; 
% name of selection table column containing label(s)
labCol = 'Tags'; 
% desired duration of clips in seconds; leave empty to maintain annotation durations, which may be variable; add 0.05s to desired clip length for Koogu training clips
winSize = 3; 
% Directory to save clips; subdirectories will be created for each unique label
clipDir = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AtlStur\ThunderClassifier_Training\KooguModel\TrainClips';
% if desired, filter clips by SNR threshold; leave empty to include all clips regardless of SNR
SNRthresh = []; 
% set to >1 to speed up audio by this factor, <1 to slow down audio
shift = 1; 

% TO DO: implement sliding window & overlap
%slidingWin = 0; % 1 if each annotation should be windowed into multiple clips (length/overlap allowing), 0 otherwise
%overlap = 0; % If using a sliding window, set the overlap in seconds

%%
% get paths and names of audio files
audioList = dir(fullfile(audioDir,['**\*',fileExt]));
audioNames = {audioList.name};
audioPaths = {audioList.folder};

% identify selection tables
selTabs = dir(fullfile(selDir,'*.txt'));

for i=1:size(selTabs,1)
    tab = readtable(fullfile(selDir,selTabs(i).name),'Delimiter',"\t",'VariableNamingRule',"preserve");
    if any(contains(tab.Properties.VariableNames,'View'))
        if (any(strmatch('Spectrogram',tab.View)) && any(strmatch('Waveform',tab.View)))
            tab = tab(strmatch('Spectrogram 1',tab.View),:);
        end
    end
    if ~any(contains(tab.Properties.VariableNames,'Delta Time (s)'))
        tab(:,'Delta Time (s)') = table(table2array(tab(:,'End Time (s)')) - table2array(tab(:,'Begin Time (s)')));
    end
    startS = table2array(tab(:,'File Offset (s)'));
    endS = startS + table2array(tab(:,'Delta Time (s)'));
    tag = table2array(tab(:,labCol));
    if isa(tag,'cell') && ~isempty(tag)
        chan = table2array(tab(:,'Channel'));
        soundFile = table2array(tab(:,'Begin File'));
        numClips = size(tab,1);

        if ~isempty(SNRthresh)
            SNRcol = find(contains(tab.Properties.VariableNames,'SNR'));
            SNR = table2array(tab(:,SNRcol));
            goodSNR = find(SNR>=SNRthresh);
            numClips = numel(goodSNR);
        else
            goodSNR = 1:numClips;
        end

        skippedClips = 0;
        savedClips = 0;

        for j=1:numClips

            soundFileInd = find(strcmp(soundFile{goodSNR(j)},audioNames));
            if ~isempty(soundFileInd)
                if contains(tag{goodSNR(j)},'/')
                    tag{goodSNR(j)} = strrep(tag{goodSNR(j)},'/','_');
                end
                if sum(isstrprop(tag{goodSNR(j)},'alphanum'))==length(tag{goodSNR(j)}) % If label string contains special characters, skip clip (directories and file names can't contain special characters)
                    info = audioinfo(fullfile(audioPaths{soundFileInd},audioNames{soundFileInd}));
                    annotDur = endS(goodSNR(j))-startS(goodSNR(j));
                    if isempty(winSize)
                        winSize = annotDur;
                    end
                    padDur = winSize - annotDur;
                    annotSampRange = round([startS(goodSNR(j)),endS(goodSNR(j))].*info.SampleRate);
                    sampRange = [annotSampRange(1)-(round(padDur*info.SampleRate/2)),annotSampRange(2)+(round(padDur*info.SampleRate/2))];

                    if (sampRange(1)<=0 || sampRange(2)>info.TotalSamples)
                        skippedClips = skippedClips + 1;
                    elseif (sampRange(2)<=0 || sampRange(2)<sampRange(1) || isempty(sampRange))
                        weirdClips = 1;
                    else
                        [waveData,Fs] = audioread(fullfile(audioPaths{soundFileInd},audioNames{soundFileInd}),sampRange,'native');
                        waveData = waveData(:,chan(goodSNR(j)));
                        if ~isempty(SNRthresh)
                            SNRstr = strrep(num2str(SNR(goodSNR(j))),'.','_');
                            saveName = [strrep(audioNames{soundFileInd},fileExt,''),'_Ch',num2str(chan(goodSNR(j))),'_',num2str(floor(startS(goodSNR(j)))),'s_',tag{goodSNR(j)},'_SNR',SNRstr,'.wav'];
                        else
                            saveName = [strrep(audioNames{soundFileInd},fileExt,''),'_Ch',num2str(chan(goodSNR(j))),'_',num2str(floor(startS(goodSNR(j)))),'s_',tag{goodSNR(j)},'.wav'];
                        end

                        if ~isfolder(fullfile(clipDir,tag{goodSNR(j)}))
                            mkdir(fullfile(clipDir,tag{goodSNR(j)}))
                        end
                        q = dir(fullfile(clipDir,tag{goodSNR(j)},'*.wav'));
                        q = {q.name};

                        if (any(strcmp(saveName,q)))
                            duplicate = 1;
                        end
                        audiowrite(fullfile(clipDir,tag{goodSNR(j)},saveName),waveData,Fs*shift,'BitsPerSample',info.BitsPerSample)
                        savedClips = savedClips+1;
                    end

                    soundFileInd = [];
                    info = [];
                    annotDur = [];
                    padDur = [];
                    annotSampRange = [];
                    sampRange = [];
                    waveData = [];
                    Fs = [];
                    saveName = [];
                else
                    fprintf('Skipping call type "%s" because of special characters\n',tag{goodSNR(j)})
                    skippedClips = skippedClips + 1;
                end
            else
                fprintf('No audio file %s\n',soundFile{goodSNR(j)})
            end
        end

        fprintf('Done with selection table %d of %d\n',i,size(selTabs,1))
        fprintf('  %d clips saved of %d selections\n',savedClips,numClips)
        if (skippedClips>0)
            fprintf('  Skipped %d selections which were too close to the start/end of their audio files, or had special characters\n',skippedClips)
        end

    elseif isa(tag,'double') & sum(isnan(tag))==length(tag)
        fprintf('Skipping selection table %d of %d because no labeled annotations present in specified column\n',i,size(selTabs,1))
    end
end