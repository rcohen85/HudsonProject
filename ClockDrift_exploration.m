%% Quantify effect of differing extents of drift on various sampling rates 
% and over different sampling periods. Assuming drift is linear (it's
% probably not, in real life).
clearvars
drift = [0.001,0.01,0.1,1]; % clock drift in s/day
Fs = [1,10,100,200]; % initial sampling rate in kHz
dur = [0,10,30,60,180,365]; % duration of deployment in days
st = datenum('2022-01-01 00:00:00');
ed = st+dur(end);

% Negative drift (clock is slowing down)
for i = 1:length(drift) % for each drift rate
    
    obsFs = zeros(length(Fs),length(dur));
    
    for j = 1:length(Fs) % for each sampling rate
        
        step = (1/(60*60*24))/(Fs(j)*1000); % desired sampling increment, in days
        nSamps = (Fs(j)*1000)*60*60*24; % number of samples per day, w no drift
        lostSamps = (Fs(j)*1000)*drift(i);% number of samples per day lost to drift
        %         % Time of first day samples if there were no drift
        %         optSpacing = st:step:(st+1-step);
        % Time of first day samples w drift
        day1Samps = linspace(st,(st+1-step),nSamps-lostSamps);
        
        %         % Amount each sample is offset from when it should've been sampled
        %         initChangePerDay = linspace(0,drift(i),(Fs(j)*60*60*24*1000));
        %         % Add offsets to get actual sample times
        %         day1Samps = optSpacing+initChangePerDay;
        lastSamp = day1Samps(end);
        
        % count samples in first second of deployment to establish starting Fs
        firstSecSamps = find(day1Samps>st & day1Samps<(st+(1/(60*60*24))));
        obsFs(j,1) = length(firstSecSamps);
        
        dy = floor(st);
        
        for k = 2:max(dur)
                thisDaySamps = linspace(lastSamp,dy+2-step,nSamps-(lostSamps*k));                
                lastSamp = thisDaySamps(end);
                dy = floor(lastSamp);
                jd = day(datetime(datestr(dy)), 'dayofyear');
                if ismember(jd,dur)
                    noonSecSamps = find(thisDaySamps>(dy+0.5) & thisDaySamps<(dy+0.5+(1/(60*60*24))-step));
                    ind = find(dur==jd);
                    obsFs(j,ind) = length(noonSecSamps);
                end
        end
        
        fprintf('Done with Fs:%dkHz\n',Fs(j));
        
    end
    
    % Plot curves
    plot(dur,obsFs,'o-')
    legend(arrayfun(@num2str,dur,'UniformOutput',0));
    
end

% Positive drift (clock is speeding up)

%% Following sections explore consequences of analyzing data at a different
% Fs ("presumed Fs") than it was actually sampled at. I assume an equal
% number of samples will always be analyzed. Therefore the time period
% these samples cover, and the true Nyquist frequency of the data, will be
% different than what would be expected based on the presumed Fs.

presumedFs = 24000;

%% Sampled at 24kHz and analyzed w correct Fs
Fs = 24000; % actual sampling rate
f = 0:((presumedFs/2))/((presumedFs/2)-1):((presumedFs/2));

% Generate signal combining tones at a few different frequencies
dt = 1/Fs; % seconds per sample
StopTime = presumedFs/Fs; % # seconds of data required to get 24k samples
t = (0:dt:StopTime)'; % seconds
F10 = sin(2*pi*10*t);
F100 = sin(2*pi*100*t);
F1000 = sin(2*pi*1000*t);
F10000 = sin(2*pi*10000*t);
combinedWave = F10+F100+F1000+F10000;

% Window & FFT
sigLength = length(combinedWave);
wind = hann(sigLength);
wWave = combinedWave.*wind;
nfft = sigLength;
waveSpec = abs(fft(wWave,nfft));
waveSpec = waveSpec(1:(sigLength-1)/2);

% Plot
figure (1)
plot(f,waveSpec)
xlim([0 presumedFs/2])
set(gca,'XScale','log');
title('FFT, True F_s: 24k Analysis F_s: 24k');

[p,freq] = pwelch(wWave,sigLength,[],[],presumedFs);
figure(2)
plot(freq,p);
set(gca,'XScale','log');
xlim([0 presumedFs/2])
title('pwelch, True F_s: 24k Analysis F_s: 24k');



%% Sampled at 23.5kHz and analyzed w wrong Fs
Fs = 23500; % actual sampling rate
f = 0:((presumedFs/2))/((presumedFs/2)-1):((presumedFs/2));

% Generate signal combining tones at a few different frequencies
dt = 1/Fs; % seconds per sample
StopTime = presumedFs/Fs; % # seconds of data required to get 24k samples
t = (0:dt:StopTime)'; % seconds
F10 = sin(2*pi*10*t);
F100 = sin(2*pi*100*t);
F1000 = sin(2*pi*1000*t);
F10000 = sin(2*pi*10000*t);
combinedWave = F10+F100+F1000+F10000;

% Window & FFT
sigLength = length(combinedWave);
wind = hann(sigLength);
wWave = combinedWave.*wind;
nfft = sigLength;
waveSpec = abs(fft(wWave,nfft));
waveSpec = waveSpec(1:(sigLength-1)/2);

% Plot
figure (3)
plot(f,waveSpec)
xlim([0 presumedFs/2])
set(gca,'XScale','log');
title('FFT, True F_s: 23.5k Analysis F_s: 24k');

[p,freq] = pwelch(wWave,sigLength,[],[],presumedFs);
figure(4)
plot(freq,p);
set(gca,'XScale','log');
xlim([0 presumedFs/2])
title('pwelch, True F_s: 23.5k Analysis F_s: 24k');

%% Sampled at 23kHz and analyzed w wrong Fs
Fs = 23000; % actual sampling rate
f = 0:((presumedFs/2))/((presumedFs/2)-1):((presumedFs/2));

% Generate signal combining tones at a few different frequencies
dt = 1/Fs; % seconds per sample
StopTime = presumedFs/Fs; % # seconds of data required to get 24k samples
t = (0:dt:StopTime)'; % seconds
F10 = sin(2*pi*10*t);
F100 = sin(2*pi*100*t);
F1000 = sin(2*pi*1000*t);
F10000 = sin(2*pi*10000*t);
combinedWave = F10+F100+F1000+F10000;

% Window & FFT
sigLength = length(combinedWave);
wind = hann(sigLength);
wWave = combinedWave.*wind;
nfft = sigLength;
waveSpec = abs(fft(wWave,nfft));
waveSpec = waveSpec(1:(sigLength-1)/2);

% Plot
figure (5)
plot(f,waveSpec)
xlim([0 presumedFs/2])
set(gca,'XScale','log');
title('FFT, True F_s: 23k Analysis F_s: 24k');

[p,freq] = pwelch(wWave,sigLength,[],[],presumedFs);
figure(6)
plot(freq,p);
set(gca,'XScale','log');
xlim([0 presumedFs/2])
title('pwelch, True F_s: 23k Analysis F_s: 24k');

%% Sampled at 25kHz and analyzed w wrong Fs
Fs = 25000; % actual sampling rate
f = 0:((presumedFs/2))/((presumedFs/2)-1):((presumedFs/2));

% Generate signal combining tones at a few different frequencies
dt = 1/Fs; % seconds per sample
StopTime = presumedFs/Fs; % # seconds of data required to get 24k samples
t = (0:dt:StopTime)'; % seconds
F10 = sin(2*pi*10*t);
F100 = sin(2*pi*100*t);
F1000 = sin(2*pi*1000*t);
F10000 = sin(2*pi*10000*t);
combinedWave = F10+F100+F1000+F10000;

% Window & FFT
sigLength = length(combinedWave);
wind = hann(sigLength);
wWave = combinedWave.*wind;
nfft = sigLength;
waveSpec = abs(fft(wWave,nfft));
waveSpec = waveSpec(1:(sigLength-1)/2);

% Plot
figure (7)
plot(f,waveSpec)
xlim([0 presumedFs/2])
set(gca,'XScale','log');
title('FFT, True F_s: 25k Analysis F_s: 24k');

[p,freq] = pwelch(wWave,sigLength,[],[],presumedFs);
figure(8)
plot(freq,p);
set(gca,'XScale','log');
xlim([0 presumedFs/2])
title('pwelch, True F_s: 25k Analysis F_s: 24k');
