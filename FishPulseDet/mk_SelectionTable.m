% Convert output from pulse detector into selection tables for viewing
% of detections in Raven. 

% directory containing original sound files
% soundDir = '';
% directory containing detector output folder(s)
baseDir = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\ShortStur\PulseDetectorOutputmetadata';
% directory where you want to save your selection table
outDir = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\ShortStur\SelectionTables\SS03\Unlabeled';
chan = 1;
% dateFormat = '_(\d{8})_(\d{6})';
% Fs = 250000; % sampling rate (Hz)
% freqRange = [8000 96000];


%%

dirSet = dir(baseDir);

for itr0 = 4:length(dirSet)
    if dirSet(itr0).isdir &&~strcmp(dirSet(itr0).name,'.')&&...
            ~strcmp(dirSet(itr0).name,'..')

        inDir = fullfile(baseDir,dirSet(itr0).name);
        thisDir = dirSet(itr0).name;
        fileSet = what(inDir);
        lfs = length(fileSet.mat);

        if ~exist(fullfile(outDir,thisDir),'dir')
            fprintf('Creating output directory %s\n',outDir)
            mkdir(fullfile(outDir,thisDir))
        end

        Selection = [];
        View = [];
        Channel = [];
        Begin_time = [];
        End_time = [];
        Low_freq = [];
        High_freq = [];
        File_offset = [];
        Begin_file = [];
        Begin_path = [];
        Begin_dateTime = [];
        End_Clocktime = [];
        Sound_type = [];
        sampleOffset = 0;

        sel = 1; % counter for selection numbers in each deployment
        fileSetStart = [];

        for itr2 = 1:lfs
            thisFile = fileSet.mat{itr2};
%             flacFile = strrep(thisFile,'.mat','.flac');
            wavFile = strrep(thisFile,'.mat','.wav');
            load(char(fullfile(inDir,thisFile)),'-mat','clickTimes','hdr',...
                'ppSignal','specClickTf','yFiltBuff','f')
            fileStart = datenum(hdr.start.dvec);
%             info = audioinfo(fullfile(soundDir,wavFile)); % audioinfo
%             gets true size of data
%             Fs = info.SampleRate;
            % reading the wave headers is less reliable, as these may be
            % incorrect in messed up files, but this is what Raven does, so
            % for annotations to show up with correct timing in Raven we
            % must be consistent
%             hdr = ioReadWavHeader(fullfile(soundDir,thisDir,wavFile),dateFormat);
%             Fs = hdr.fs;
            
%             if itr2==1
%                 fileSetStart = fileStart;
%             end

            if ~isempty(clickTimes)

                [~,keepers] = unique(clickTimes(:,1));
                clickTimes = clickTimes(keepers,:);
                numDets = size(clickTimes,1);
%                 posDnum = (clickTimes(:,1)/(60*60*24)) + fileStart;
%                 posDvec = datevec(datetime(posDnum,'ConvertFrom','datenum'));
%                 dStrings = cellstr([num2str(posDvec(:,1)),repmat('/',numDets,1),...
%                     num2str(posDvec(:,2)),repmat('/',numDets,1),...
%                     num2str(posDvec(:,3)),repmat(' ',numDets,1),...
%                     num2str(posDvec(:,4)),repmat(':',numDets,1),...
%                     num2str(posDvec(:,5)),repmat(':',numDets,1),num2str(posDvec(:,6),'%-.9f')]);
%                 clockEnd = cellstr(strcat(num2str(posDvec(:,4)),repmat(':',numDets,1),...
%                     num2str(posDvec(:,5)),repmat(':',numDets,1),num2str(posDvec(:,6),'%-.9f')));

%                 fileSet_offset = (fileStart - fileSetStart)*(60*60*24); % seconds since beginning of file set
                

%                 % Calculate -10dB (-90% power) bandwidth edges for annotation bounding boxes
%                 fLow = [];
%                 fHigh = [];
%                 for itr3 = 1:length(clickTimes)
%                     [valMx, posMx] = max(specClickTf(itr3,:)); % find amplitude and value of peak freq
%                     low = valMx-10; % -10dB value relative to peak
%                     smoothSpec = fastsmooth(specClickTf(itr3,:),10);
% 
%                     % from peak freq, walk along spectrogram until low is reached on either side
%                     slopeup=fliplr(smoothSpec(1,1:posMx));
%                     for e10dB=1:length(slopeup)
%                         if slopeup(e10dB)<low %stop at value < -3dB: point of lowest frequency
%                             break
%                         end
%                     end
%                     slopedown=smoothSpec(1,posMx:length(specClickTf(itr3,:)));
%                     for o10dB=1:length(slopedown)
%                         if slopedown(o10dB)<low %stop at value < -3dB: point of highest frequency
%                             break
%                         end
%                     end
% 
%                     %calculation from spectrogram -> from 0 to 100kHz in 256 steps (FFT=512)
%                     high10dB = f(posMx+o10dB-1)*1000; %(Fs/(2))*((posMx+o10dB)/length(specClickTf(itr3,:))); %-10dB highest frequency in Hz
%                     low10dB = f(posMx-e10dB+1)*1000; %(Fs/(2))*((posMx-e10dB)/length(specClickTf(itr3,:))); %-10dB lowest frequency in Hz
% 
%                     fLow = [fLow;low10dB];
%                     fHigh = [fHigh;high10dB];

                % Calculate annotation bounding box frequency limits
                % (highest/lowest freqs within freqRange where the power attains 10% of the peak amplitude)
                fLow = [];
                fHigh = [];
                for itr3 = 1:size(clickTimes,1)
                    [valMx, ~] = max(specClickTf(itr3,:)); % find amplitude of peak freq
                    low = valMx-10; % -10dB value relative to peak
                    
                    highInd = find(specClickTf(itr3,:)>=low);

                    low10dB = f(highInd(1))*1000;
                    high10dB = f(highInd(end))*1000;

                    fLow = [fLow;low10dB];
                    fHigh = [fHigh;high10dB];
                end

                Selection = [Selection;(sel:(sel+numDets-1))'];
                View = [View;cellstr(repmat('Waveform 1',numDets,1))];
                Channel = [Channel;repmat(chan,numDets,1)];
%                 Begin_time = [Begin_time;clickTimes(:,1)+fileSet_offset];  % start position of pulses relative to beginning of file set (s)
%                 End_time = [End_time;clickTimes(:,2)+fileSet_offset]; % end position of pulses relative to beginning of file set (s)
%                 Begin_time = [Begin_time;((clickTimes(:,1)*Fs)+sampleOffset)./Fs];  % start position of pulses relative to beginning of file set (s)
%                 End_time = [End_time;((clickTimes(:,2)*Fs)+sampleOffset)./Fs]; % end position of pulses relative to beginning of file set (s)
                Begin_time = [Begin_time;clickTimes(:,1)];  % start position of pulses relative to beginning of file (s)
                End_time = [End_time;clickTimes(:,2)]; % end position of pulses relative to beginning of file (s)
                Low_freq = [Low_freq;fLow]; % lower freq bound for annotation box
                High_freq = [High_freq;fHigh]; % upper freq bound for annotation box
                File_offset = [File_offset;clickTimes(:,1)]; % position of pulses relative to beginning of file (s)
                Begin_file = [Begin_file;repmat(wavFile,numDets,1)]; % sound file each pulse is found in
%                 Begin_path = [Begin_path;repmat(soundDir,numDets,1)]; % path to where sound files are found
%                 Begin_dateTime = [Begin_dateTime;dStrings]; % absolute timestamps of pulses
%                 End_Clocktime = [End_Clocktime;clockEnd];  % pulse end clock time (no calendar year/month/day)
%                 Sound_type = [Sound_type;repmat('Pulse',numDets,1)];

                sel = sel+numDets;
            end

            fprintf('Done with file %d of %d \n',itr2,lfs)

%             if (size(Begin_time,1)>= 10000 || itr2 == lfs)

%                 saveTable = table(Selection,View,Channel,Begin_time,End_time,Low_freq,...
%                     High_freq,File_offset,Begin_file,Begin_path,Begin_dateTime,End_Clocktime,Sound_type,...
%                     'VariableNames',{'Selection','View','Channel','Begin Time (s)',...
%                     'End Time (s)','Low Freq (Hz)','High Freq (Hz)','File Offset (s)','Begin File',...
%                     'Begin Path','Begin Date Time','End Clock Time','Sound Type'});

%                 saveTable = table(Selection,View,Channel,Begin_time,End_time,Low_freq,...
%                     High_freq,File_offset,Begin_file,Begin_path,Sound_type,...
%                     'VariableNames',{'Selection','View','Channel','Begin Time (s)',...
%                     'End Time (s)','Low Freq (Hz)','High Freq (Hz)','File Offset (s)','Begin File',...
%                     'Begin Path','Sound Type'});

                saveTable = table(Selection,View,Channel,Begin_time,End_time,Low_freq,...
                    High_freq,Begin_file,File_offset,...
                    'VariableNames',{'Selection','View','Channel','Begin Time (s)',...
                    'End Time (s)','Low Freq (Hz)','High Freq (Hz)','Begin File','File Offset (s)'});

                if itr2 == lfs 
                    fprintf('Done with directory %d of %d \n',itr0,length(dirSet))

                end
                
%                 firstDet = datestr(Begin_dateTime{1,1},'yyyymmdd_HHMMSS');
%                 lastDet = datestr(Begin_dateTime{end,1},'yyyymmdd_HHMMSS');
%                 saveName =  [fullfile(outDir,dirSet(itr0).name),'_SelectionTable_',firstDet,'_',lastDet,'.txt'];
                saveName = fullfile(outDir,thisDir,strrep(thisFile,'.mat','.txt'));
                writetable(saveTable,saveName,'Delimiter','\t');

                Selection = [];
                View = [];
                Channel = [];
                Begin_time = [];
                End_time = [];
                Low_freq = [];
                High_freq = [];
                File_offset = [];
                Begin_file = [];
                Begin_path =[];
                Begin_dateTime = [];
                End_Clocktime = [];
                Sound_type = [];
                firstDet = [];
                lastDet = [];
%                 sel = 1;
                

%             end
            
%             sampleOffset = sampleOffset + info.TotalSamples;
%               sampleOffset = sampleOffset + (hdr.Chunks{hdr.dataChunk}.DataSize/...
%                   (hdr.Chunks{hdr.fmtChunk}.Info.nBytesPerSample * hdr.Chunks{hdr.fmtChunk}.Info.nChannels));

        end

    end
end