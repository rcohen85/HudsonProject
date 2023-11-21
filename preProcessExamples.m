inFile = 'C:\Users\rec297\CCB\Teaching\NYSDEC01_20171101_151140.wav';
info = audioinfo(inFile);
[y,Fs] = audioread(inFile);
y = double(y(:,1));

NFFT = 256;
win = hann(NFFT);
overlap = floor(NFFT*0.75);

[s,f,t] = spectrogram(y(1:60*info.SampleRate),win,overlap,NFFT,info.SampleRate);
s = 20*log10(abs(s));

figure(2),clf
subplot(2,1,1)
plot(y(1:60*info.SampleRate))
xlabel('Sample')
ylabel('Normalized Amplitude')
xlim([1,length(y(30*info.SampleRate:90*info.SampleRate))])
ylim([-0.05 0.05])
set(gca,'FontSize',14)
subplot(2,1,2)
colormap jet
imagesc(t,f,s)
set(gca,'YDir','normal')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
set(gca,'FontSize',14)
clim([-60,10])

saveas(gcf,'C:\Users\rec297\CCB\Teaching\Fin20HzPulses_lowSNR.png')


%% Bandpass filter timeseries, denoise spectrogram

[b,a] = butter(5,[0.0375 0.0875]);
cleanData = filter(b,a,y);

M = movmean(s,60,2);
s_norm = s-M;
s_norm(s_norm<0) = 0;

figure(2),clf
subplot(2,1,1)
plot(cleanData(1:60*info.SampleRate))
xlabel('Sample')
ylabel('Normalized Amplitude')
xlim([1,60*info.SampleRate])
ylim([-0.007 0.007])
set(gca,'FontSize',14)
subplot(2,1,2)
colormap jet
imagesc(t,f,s_norm)
set(gca,'YDir','normal')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
set(gca,'FontSize',14)
clim([7 22])

saveas(gcf,'C:\Users\rec297\CCB\Teaching\preProcessedData.png')