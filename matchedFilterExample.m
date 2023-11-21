% Create synthetic call waveform template
fs = 400; % Sampling frequency
dt = 1/fs; % seconds per sample 
dur = 2; % duration of signal (seconds)
t = (0:dt:dur)'; % timestamps of samples (seconds)
F = 20; % desired signal frequency (Hz) 
amp = [0.*(t(t<0.5));(t(t>=0.5 & t<0.9)-0.5).*0.2;0.08.*(ones(sum(t>=0.9 & t<=1.1),1));(flipud(t(t>=1.1 & t<1.5))-1.1).*0.2;0.*(t(t>1.5))];
data = amp.*sin(2*pi*F*t);

figure(1);clf
plot(data,'LineWidth',1.5);
xlim([1,length(data)])
ylim([-0.1 0.1])
xlabel('Sample')
ylabel('Amplitude')
set(gca,'FontSize',14)
title('Synthetic 20Hz Pulse Template');

saveas(gcf,'C:\Users\rec297\CCB\Teaching\MatchedFilterSyntheticTemplate.png')

figure(1);clf
plot(data,'-r');
xlim([1,length(data)])
ylim([-0.1 0.1])
xticklabels({})
yticklabels({})
set(gca,'color','none')
export_fig('C:\Users\rec297\CCB\Teaching\MatchedFilterSyntheticTemplate_noBackground.png','-transparent',1)

%% Load real data
inFile = 'C:\Users\rec297\CCB\Teaching\NYSDEC01_20171101_151140.wav';
info = audioinfo(inFile);
[y,Fs] = audioread(inFile);
y = double(y(:,3));

% Get call to use as natural waveform template
niceCall = y(37.3*info.SampleRate:39.3*info.SampleRate);
[b,a] = butter(5,0.15);
cleanCall = filter(b,a,niceCall);

figure(1);clf
plot(cleanCall,'LineWidth',1.5);
ylim([-0.1 0.1])
xlim([1,length(cleanCall)])
xlabel('Sample')
ylabel('Amplitude')
set(gca,'FontSize',14)
title('Natural 20Hz Pulse Template');

saveas(gcf,'C:\Users\rec297\CCB\Teaching\MatchedFilterNaturalTemplate.png')

%% Plot raw data w 20Hz pulses

NFFT = 256;
win = hann(NFFT);
overlap = floor(NFFT*0.75);

[s,f,t] = spectrogram(y(30*info.SampleRate:90*info.SampleRate),win,overlap,NFFT,info.SampleRate);
s = 20*log10(abs(s));

figure(2),clf
subplot(2,1,1)
plot(y(30*info.SampleRate:90*info.SampleRate))
xlabel('Sample')
ylabel('Normalized Amplitude')
xlim([1,length(y(30*info.SampleRate:90*info.SampleRate))])
ylim([-0.1 0.12])
set(gca,'FontSize',14)
subplot(2,1,2)
colormap jet
imagesc(t,f,s)
set(gca,'YDir','normal')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
set(gca,'FontSize',14)
clim([-35,15])

saveas(gcf,'C:\Users\rec297\CCB\Teaching\Fin20HzPulses.png')


%% Cross-correlate waveform template w raw waveform data

% % fancy way (same results as simple way)
% synthFilt = data(end:-1:1);
% natFilt = cleanCall(end:-1:1);
% y_synthFilt = filter(synthFilt,1,y(30*info.SampleRate:90*info.SampleRate));
% y_natFilt = filter(natFilt,1,y(30*info.SampleRate:90*info.SampleRate));
% 
% figure(1),clf
% plot(y_natFilt)
% xlabel('Lag')
% ylabel('Cross-Correlation')
% set(gca,'FontSize',14)

% simple way
[c,lags] = xcorr(y(30*info.SampleRate:90*info.SampleRate),data);


figure(5),clf
plot(lags,c)
yline(0.2,'-r','LineWidth',1.5)
xlabel('Lag')
ylabel('Cross-Correlation')
xlim([0 lags(end)])
set(gca,'FontSize',14)

saveas(gcf,'C:\Users\rec297\CCB\Teaching\MatchedFilterOutput.png')


%%  Zoomed spectrogram showing downsweep contour

inFile = 'C:\Users\rec297\CCB\Teaching\NYSDEC01_20171101_151140.wav';
info = audioinfo(inFile);
[y,Fs] = audioread(inFile,[30*info.SampleRate,90*info.SampleRate]);
y = double(y(:,3));

NFFT = 256;
win = hann(NFFT);
overlap = floor(NFFT*0.75);

[s,f,t] = spectrogram(y,win,overlap,NFFT,info.SampleRate);
s = 20*log10(abs(s));

figure(5),clf
colormap jet
imagesc(t,f,s)
set(gca,'YDir','normal')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
set(gca,'FontSize',14)
clim([-35,15])
ylim([0 50])
% xlim([4 31])

saveas(gcf,'C:\Users\rec297\CCB\Teaching\Fin20HzPulseSpecZoom.png')


%% Spectrogram template

inFile = 'C:\Users\rec297\CCB\Teaching\NYSDEC01_20171101_151140.wav';
info = audioinfo(inFile);
[y,Fs] = audioread(inFile);
y = double(y(:,3));

% Get call to use as specctrogram template
niceCall = [y(36*info.SampleRate:41*info.SampleRate)];
[b,a] = butter(5,0.15);
cleanCall = filter(b,a,niceCall);

NFFT = 256;
win = hann(NFFT);
overlap = floor(NFFT*0.75);

[s,f,t] = spectrogram(cleanCall,win,overlap,NFFT,info.SampleRate);
s = 20*log10(abs(s));
s(s<-25) = NaN;

figure(1)
colormap jet
imagesc(t,f,s)
set(gca,'YDir','normal')
xticklabels({})
yticklabels({})
ylim([0 50])
clim([-25,10])

saveas(gcf,'C:\Users\rec297\CCB\Teaching\SpectrogramTemplate.png')

%% Perform spectrogram cross-correlation

inFile = 'C:\Users\rec297\CCB\Teaching\NYSDEC01_20171101_151140.wav';
info = audioinfo(inFile);
[y,Fs] = audioread(inFile);
y = double(y(:,1));

% Get call to use as natural spectrogram template
niceCall = y(37.3*info.SampleRate:39.3*info.SampleRate);
[b,a] = butter(5,0.15);
cleanCall = filter(b,a,niceCall);

y = y(30*info.SampleRate:90*info.SampleRate);

win = 801;
NFFT = 256;
noverlap = floor(win*0.75);

S = [];
T = [];
k = 1;
for i=1:(length(y)-win)
% s = xspectrogram(cleanCall,y(i:(i+win-1)),win,noverlap,NFFT);
s = xspectrogram(cleanCall,cleanData(i:(i+win-1)),win,noverlap,NFFT);

S = [S;max(s(2:end))];
end

figure(5),clf
plot(S,'-','LineWidth',1.5)
yline(0.2e-4,'-r','LineWidth',1.5)
xlabel('Lag')
ylabel('Cross-Correlation')
xlim([0 length(S)])
set(gca,'FontSize',14)

saveas(gcf,'C:\Users\rec297\CCB\Teaching\SpecXCorrOutput.png')


