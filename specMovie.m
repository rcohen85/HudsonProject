
inFile = 'C:\Users\rec297\Downloads\20230603_array2_cleaned.wav';
outFile = 'C:\Users\rec297\CCB\MDWaveVid_Hotter.avi';
frameLen = 10; % duration of each spectrogram in s
NFFT = 4000;

info = audioinfo(inFile);
nSamps = frameLen*info.SampleRate;
frames = ceil(info.TotalSamples/nSamps);
win = hann(NFFT);
noverlap = floor(NFFT*0.75);

vid = VideoWriter(outFile,'Motion JPEG AVI');
vid.FrameRate = 2;
open(vid);
ind=1;

for i=1:frames
    if i<frames
    x = double(audioread(inFile,[ind,(ind+nSamps-1)],'native'));
    elseif i==frames
        x = double(audioread(inFile,[ind,info.TotalSamples],'native'));
    end
    figure(1)
    spectrogram(x,win,noverlap,NFFT,info.SampleRate,'yaxis')
    ylim([0,3])
    colormap jet
    clim([-50,35])
    colorbar('off')
    F(i) = getframe(gcf);
    writeVideo(vid,F(i));
    ind = ind+nSamps;
end
close(vid);
toc