%% For a few irregularly-spaced playbacks, specify start time of each playback

% inDir = '/Users/rec297/Documents/CCB/UnderwaterPlaybackTest';
% outDir = '/Users/rec297/Documents/CCB/UnderwaterPlaybackTest';
% BD = 16; % bit depth
% toneStart = [1400.4,1402.7,1404.75,1406.85,1408.92,1411.02,1413.02];
% toneDur = 0.2; % s
% clipDur = 0.1;
% freqs = [150,300,1000,2000,3000];
% refPres = 1e-6; % reference pressure for SPL measurements (1 muPa for underwater, 20 muPa for in air)
% cal = 176.6;
% cal_lin = power(10, cal/20);
%
% fileList = dir(fullfile(inDir,'*.wav'));
%
% for k=1:length(fileList)
%     soundFile = fileList(k).name;
%     fileStart = datetime(regexp(soundFile,'.\d{14}','match'),'InputFormat','.uuuuMMddHHmmss');
%     [y,Fs] = audioread(fullfile(inDir,soundFile));
%     y = double(y);
%     leftover = toneDur-clipDur;
%
%     startTime = [];
%     endTime = [];
%
%     tones = [];
%     testSPLs = [];
%     for j=1:size(toneStart,2)
%
%         counter = 0;
%
%         for i=1:(length(freqs))
%
%             clipStart = toneStart(k,j)+(leftover/2)+counter;
%             clipEnd = clipStart+clipDur;
%
%             clip = y(floor(clipStart*Fs):floor(clipEnd*Fs));
%             if ~isempty(tones)
%                 if length(clip)>size(tones,2)
%                     clip = clip(1:size(tones,2));
%                 elseif length(clip)<size(tones,2)
%                     clip = [clip;clip(end)];
%                 end
%             end
%             tones = [tones;clip'];
%             startTime = [startTime;clipStart];
%             endTime = [endTime;clipEnd];
%             counter = counter+toneDur;
%
%         end
%
%     end
%
%     % Make selection table to verify data are being pulled from appropriate times
%     Selection = 1:length(startTime);
%     View = repmat('Spectrogram',length(startTime),1);
%     Channel = ones(length(startTime),1);
%     lowFreq = zeros(length(startTime),1);
%     highFreq = (Fs/2)*ones(length(startTime),1);
%
%     [path,name,ext] = fileparts(soundFile);
%     beginPath = repmat(fullfile(inDir,[name,ext]),length(startTime),1);
%     beginsoundFile = repmat([name,ext],length(startTime),1);
%
%     saveTable = table(Selection',View,Channel,startTime,endTime,lowFreq,...
%         highFreq,beginPath,beginsoundFile,startTime,...
%         'VariableNames',{'Selection','View','Channel','Begin Time (s)',...
%         'End Time (s)','Low Freq (Hz)','High Freq (Hz)','Begin Path','Begin File','File Offset (s)'});
%     saveName = fullfile(outDir,strrep(soundFile,ext,'.txt'));
%     writetable(saveTable,saveName,'Delimiter','\t');
%
% end
%
%
% tones = tones-mean(tones,2); % remove DC offset
% nfft = size(tones,2);
% specs = abs(fft(tones',nfft))';
%
% % get back to units of ADC counts
% specs = specs./(nfft/2);
%
% % reduce to positive frequencies
% specs = specs(:,1:floor(nfft/2)+1);
%
% f = Fs/nfft*(0:(nfft/2));
% % plot(f,specs(4,:)), ylim([1,((2^16)/2)-1])
%
% %     [amps,mInd] = max(specs,[],2); % don't do this, sometimes text tones were not max amplitudes
% %     maxFreqs = f(mInd);
% mInd = [];
% freqVec = repmat(freqs,1,3);
% for l=1:length(freqVec)
%     [M,ind] = min(abs(freqVec(l)-f));
%     mInd = [mInd;ind];
% end
% % mInd = repelem(mInd,3);
% % missMeasurements = length(mInd)-size(tones,1);
% % mInd = mInd(1:(end-missMeasurements)); % Sound level meter stopped logging early, don't have SPL measurements for last few tones
% testFreqs = f(mInd);
% allData(k,1).TestFreqs = testFreqs;
% amps = specs(sub2ind([size(specs,1),size(specs,2)],(1:size(specs,1))',mInd)); % why is this so complicated?
%
% % Look at received level at each tested frequency
% RL_dB = 20*log10((amps.*cal_lin)./refPres);
%
% figure(2),clf
% ind=1;
% hold on
% for i=1:3
%     plot(freqs,RL_dB(ind:ind+length(freqs)-1),'o-')
%     ind=ind+length(freqs)
% end
% hold off
% legend({'Full Volume','Half Volume','Lowest Volume'})
% set(gca,'FontSize',14)


%% For regularly-spaced playbacks, perhaps across several sound files,
% specify start time of first playback in each sound file

inDir = 'W:\projects\2024_NYSDEC_NY_093213\Sounds\20240614_Playbacks_TowedST600\FilesWTones';
outDir = 'W:\projects\2024_NYSDEC_NY_093213\20240614_PlaybackSelectionTables';
BD = 16; % bit depth
toneStart = [1438.97,0.836,1.12,0.4288];
toneDur = 0.2; % s
clipDur = 0.1; % s
freqs = [150,300,1000,2000,3000];
gap = 1; % s between playback iterations
refPres = 1e-6; % reference pressure for SPL measurements (1 muPa for underwater, 20 muPa for in air)
cal = 176.6;
cal_lin = power(10, cal / 20);


fileList = dir(fullfile(inDir,'*.wav'));
toneStarts = [];
toneEnds = [];
tones = [];
for k=4:length(fileList)
    soundFile = fileList(k).name;
    fileStart = datetime(regexp(soundFile,'.\d{12}','match'),'InputFormat','.uuuuMMddHHmmss');
    [y,Fs] = audioread(fullfile(inDir,soundFile));
    info = audioinfo(fullfile(inDir,soundFile));
    y = double(y);

    leftover = toneDur-clipDur;

    startTime = [];
    endTime = [];
    clips = [];
    clipStart = 0;
    clipEnd = 0;
    counter = 0;

    while floor(clipEnd)<=info.TotalSamples
        for i=1:length(freqs)
            clipStart = floor((toneStart(k)+(leftover/2)+counter)*Fs);
            clipEnd = clipStart+(Fs*clipDur)-1;
            if k==2 & 515<clipEnd/Fs & clipEnd/Fs<1438
                clipStart = clipStart + floor(9.84*Fs);
                clipEnd = clipEnd + floor(9.84*Fs);
            elseif k==2 & clipEnd/Fs>1438
                clipStart = clipStart + floor(10.29*Fs);
                clipEnd = clipEnd + floor(10.29*Fs);
            elseif k==3 && clipEnd/Fs>214.5
                clipStart = clipStart + 58.95*Fs;
                clipEnd = clipEnd + 58.95*Fs;
            elseif k==4 && clipEnd/Fs>526
                clipStart = clipStart + 23.43*Fs;
                clipEnd = clipEnd + 23.43*Fs;
            end
            if floor(clipEnd)<=info.TotalSamples
                clip = y(clipStart:clipEnd);
                clips = [clips;clip'];
                startTime = [startTime;clipStart/Fs];
                endTime = [endTime;clipEnd/Fs];
                counter = counter+toneDur;
            end
        end
        counter = counter+gap;
    end

    tone = [tones;clips];
    toneStarts = [toneStarts;startTime];
    toneEnds = [toneEnds;endTime];


    % Make selection table to verify data are being pulled from appropriate times
    Selection = 1:length(startTime);
    View = repmat('Spectrogram',length(startTime),1);
    Channel = ones(length(startTime),1);
    lowFreq = zeros(length(startTime),1);
    highFreq = (Fs/2)*ones(length(startTime),1);

    [path,name,ext] = fileparts(soundFile);
    beginPath = repmat(fullfile(inDir,[name,ext]),length(startTime),1);
    beginsoundFile = repmat([name,ext],length(startTime),1);

    saveTable = table(Selection',View,Channel,startTime,endTime,lowFreq,...
        highFreq,beginPath,beginsoundFile,startTime,...
        'VariableNames',{'Selection','View','Channel','Begin Time (s)',...
        'End Time (s)','Low Freq (Hz)','High Freq (Hz)','Begin Path','Begin File','File Offset (s)'});
    saveName = fullfile(outDir,strrep(soundFile,ext,'.txt'));
    writetable(saveTable,saveName,'Delimiter','\t');


end

% Find amplitude of each clip at test tone freq
tones = tones-mean(tones,2); % remove DC offset
nfft = size(tones,2);
specs = abs(fft(tones',nfft))';

% get back to units of ADC counts
specs = specs./(nfft/2);

% reduce to positive frequencies
specs = specs(:,1:(nfft/2)+1);

f = Fs/nfft*(0:(nfft/2));

mInd = [];
for l=1:length(freqVec)
    [M,ind] = min(abs(freqVec(l)-f));
    mInd = [mInd;ind];
end
mInd = repelem(mInd,3);
missMeasurements = length(mInd)-size(tones,1);
mInd = mInd(1:(end-missMeasurements)); % Sound level meter stopped logging early, don't have SPL measurements for last few tones
testFreqs = f(mInd);
allData(k,1).TestFreqs = testFreqs;
amps = specs(sub2ind([size(specs,1),size(specs,2)],(1:size(specs,1))',mInd)); % why is this so complicated?

% Look at response at each tested frequency

testPa = ((10.^(testSPLs./20))*refPres); % Convert received SPLs to Pa
sysResponse = amps./testPa; % ADC counts per Pa
norm_sysResponse = sysResponse/(((2^BD)/2)-1); % percent of full scale per Pa
fullScale_Pa = 1./norm_sysResponse; % Pa resulting in full scale
fullScale_dB = 20*log10(fullScale_Pa/refPres); % dB re refPres resulting in full scale

allData(k,1).Response = norm_sysResponse;
allData(k,1).FullScale_dB = fullScale_dB;
allData(k,1).Playback_dB = testSPLs;

g = findgroups(testFreqs);
%     meanTestFreqs = splitapply(@mean,testFreqs,g);
%     allData(k,1).TestFreqs = testFreqs;

meds1 = splitapply(@median,norm_sysResponse,g');
%     allData(k,1).Response = meds1;

figure(1),clf
boxplot(20*log10(norm_sysResponse),testFreqs) % converting system response to dB FS (20*log10(sysResponse/maxADCcounts))
h = gca;
hold on
plot(h.XTick,20*log10(meds1),'-')
xlabel('Frequency (Hz)')
ylabel('dBFS @ 1 Pa / 94 dB')
ylim([-120 0])
title([dev,' Frequency Response'])
set(gca,'FontSize',13)
hold off

saveName = fullfile(outDir,[dev,'_FrequencyResponse.png']);
saveas(gcf,saveName,'png');

meds2 = splitapply(@median,fullScale_dB,g');
%     allData(k,1).FullScale_dB = meds2;

figure(2),clf
boxplot(fullScale_dB,testFreqs)
h = gca;
hold on
plot(h.XTick,meds2,'-')
xlabel('Frequency (Hz)')
ylabel('Amplitude (dB re 20 \muPa)')
ylim([90 220])
title([dev,' Received Levels Eliciting Full Scale Response'])
set(gca,'FontSize',13)
hold off

saveName = fullfile(outDir,[dev,'_FullScaleRL.png']);
saveas(gcf,saveName,'png');

meds3 = splitapply(@median,testSPLs,g');
%     allData(k,1).Playback_dB = meds3;

figure(3),clf
boxplot(testSPLs,testFreqs)
h = gca;
hold on
plot(h.XTick,meds3,'-')
xlabel('Frequency (Hz)')
ylabel('Amplitude (dB re 20 \muPa)')
title([dev,' Playback Received Levels'])
set(gca,'FontSize',13)
hold off

saveName = fullfile(outDir,[dev,'_PlaybackRL.png']);
saveas(gcf,saveName,'png');
