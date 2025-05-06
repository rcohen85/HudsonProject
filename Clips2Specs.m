% Generate spectrograms of a set of clips with user-defined spectrogram
% parameters. Save spectrogram values as .txt files and also images as
% .png files

% Directory containing clips
inDir = 'P:\users\cohen_rebecca_rec297\CCB\DCLDE2024\MA_Data\Train\Clips\3s_Centered_StaticWindow\NARW';
% File extension of clips
fileExt = '.wav';
% Directory to save spectrograms (data and images)
saveDir = 'P:\users\cohen_rebecca_rec297\CCB\Test';
% Spectrogram settings
nfft = 350; % FFT length
overlap = 0.85; % percent overlap between data segments
minFreq = 50; % minimum frequency (Hz)
maxFreq = 300; % maximum frequency (Hz)

%%

clipList = dir(fullfile(inDir,['*',fileExt]));

if ~isfolder(saveDir)
    mkdir(saveDir)
end

for i=1:length(clipList)
    
    % Read in clip
    [y, Fs] = audioread(fullfile(inDir,clipList(i).name));

    % Create window function
    win = hann(nfft);

    % Compute spectrogram
    [s,f,t] = spectrogram(y,win,round(nfft*overlap),nfft,Fs);
    s = 20*log10((abs(s)*2).^2);

    % Trim spectrogram to desired frequency range
    [m,minInd] = min(abs(f-minFreq));
    [m,maxInd] = min(abs(f-maxFreq));
    s = s(minInd:maxInd,:);

    % Plot and save image of spectrogram
    figure(1)
    colormap jet
    imagesc(t,f(minInd:maxInd),s)
    set(gca,'YDir','normal')
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')

    % Save spectrogram values
    specTab = table(s);
    specName = fullfile(saveDir,strrep(clipList(i).name,fileExt,'.txt'));
    writetable(specTab,specName,'WriteVariableNames',false,'Delimiter','\t');

end
