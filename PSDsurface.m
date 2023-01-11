% Take PSD estimates from Triton Soundscape analysis and plot as 2D and 3D surfaces to see changes in soundscape over time

psd = readtable('W:\projects\2022_NOAA-NERRS_HudsonNY_144488\SoundscapeMetrics\BlackCreek\BC05\BC05_PSD_1h.csv');

timestamps = datenum(datetime(table2array(psd(:,1)),'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z'));

% get frequency and time vectors, and PSD values, as double
freqs = psd.Properties.VariableNames(2:end);
freqs = str2double(cellfun(@(x) extractAfter(x,'PSD_'),freqs,'UniformOutput',0));
specs = table2array(psd(:,2:end));

figure(99)
imagesc(timestamps,freqs,specs')
datetick('x')
xlim([min(timestamps),max(timestamps)])
set(gca,'YDir','normal')
xlabel('Date')
ylabel('Frequency (Hz)')

figure(999)
surf(timestamps,freqs,specs','EdgeColor','none')
datetick('x')
xlim([min(timestamps),max(timestamps)])
xlabel('Date')
ylabel('Frequency (Hz)')


% if necessary, reduce size of PSD array by time averaging and
% downsampling frequency resolution
%
% average n-sized time chunks
% n = 5;
% % maxN = (size(specs,1)-rem(size(specs,1),n))-floor(rem(size(specs,1),n)/2); % might not always want this "+1"?
% smoothSpecs = [];
% smoothTimestamps = [];
% i=1;
% j=1;
% while i <= size(specs,1)-(n-1)
% smoothSpecs(j,:) = mean(specs(i:(i+n-1),:),1);
% smoothTimeStamps(j) = timestamps(i+((n-1)/2));
% i=i+n;
% j=j+1;
% end
% leftover = size(specs,1)-i+1;
% smoothSpecs(j,:) = mean(specs(i:end,:),1);
% smoothTimeStamps(j) = timestamps(i+floor(leftover/2));
% 
% % downsample frequency resolution
% dsSpecs = downsample(smoothSpecs',4)';
% dsFreqs = downsample(freqs,4);
% 
% figure(1)
% imagesc(smoothTimeStamps,dsFreqs,dsSpecs')
% datetick('x')
% xlim([min(smoothTimeStamps),max(smoothTimeStamps)])
% set(gca,'YDir','normal')
% xlabel('Date')
% ylabel('Frequency (Hz)')
% 
% figure(2)
% surf(smoothTimeStamps,dsFreqs,dsSpecs','EdgeColor','none')
% datetick('x')
% xlim([min(smoothTimeStamps),max(smoothTimeStamps)])
% xlabel('Date')
% ylabel('Frequency (Hz)')


