function [ freqTable ] = getBandTable(fftBinSize, bin1CenterFrequency, fs, base, ...
    bandsPerDivision, firstOutputBandCenterFrequency, useFFTResAtBottom)
% getBandTable(): This is generic software that returns a three column array:
%     with the start, center, stop frequecies for logarthimically spaced
%     frequency bands such as millidecades or decidecades(third octaves base %     10) or third octaves base 2. These tables are passed to
%     'getBandSquaredSoundPressure' to convert square spectra to band
%     levels. The squared pressure can be converted to power spectral
%     density by dividing by the bandwidths.
% Inputs:
%       fftBinSize - the size of the FFT bins in Hz that subsequent
%       	processing of data will use.
%       bin1CenterFrequency - this is the center frequency in Hz of the FFT
%           spectra that will be passed to subsequent processing, normally
% 		this should be zero.
%       fs - data sampling frequency in Hz
%       base - base for the band levels, generally 10 or 2.
%       bandsPerDivision - the number of bands to divide the spectrum into
%           per increase by a factor of 'base'. A base of 2 and
%           bandsPerDivision of 3 results in third octaves base 2. Base 10
%           and bandsPerDivision of 1000 results in milliDecades.
%       firstOutputBandCenterFrequency: this is the frquency where the output
%           bands will start.
%       useFFTResAtBottom: In some cases, like milliDecades, we do not want
%           to have logarithmically spaced frequency bands across the full
%           spectrum, instead we have the option to have bands that are equal
%           FFTBinSize. The switch to log spacing is made at the band that
%           has a bandwidth greater than FFTBinSize and such that the
%           frequency space between band center frequencies is at least
%            FFTBinSize.
% Outputs:
%       Three column array where column 1 is the lowest frequency of the
%       band, column 2 is the center frequency, and 3 is the highest
%       frequency
% Author: Bruce Martin, JASCO Applied Sciences, Feb 2020; updated April 2021
%          bruce.martin@jasco.com.

bandCount = 0;

maxFreq = fs/2;

lowSideMultiplier = power(base, -1/(2*bandsPerDivision));
highSideMultiplier = power(base, 1/(2*bandsPerDivision));

% count the number of bands:
logBinCount = 0;
centerFreq = 0;
if (useFFTResAtBottom)
    binWidth = 0;
    while (binWidth < fftBinSize)
        bandCount = bandCount + 1;
        centerFreq =  getCenterFreq(base, bandsPerDivision, bandCount,                firstOutputBandCenterFrequency);
        binWidth = highSideMultiplier*centerFreq - lowSideMultiplier*centerFreq;
    end
    % now keep counting until the difference between the log spaced
    % center frequency and new frequency is greater than .025
    centerFreq =  getCenterFreq(base, bandsPerDivision, bandCount, firstOutputBandCenterFrequency);
    linearBinCount = round(centerFreq / fftBinSize);
    dC = abs(linearBinCount * fftBinSize - centerFreq) + .1;
    while (abs(linearBinCount * fftBinSize - centerFreq) < dC)
        dC = abs(linearBinCount * fftBinSize - centerFreq);
        bandCount = bandCount + 1;
        linearBinCount = linearBinCount + 1;
        centerFreq =  getCenterFreq(base, bandsPerDivision, bandCount,firstOutputBandCenterFrequency);
    end
    linearBinCount = linearBinCount - 1;
    bandCount = bandCount - 1;

    if (fftBinSize * linearBinCount > maxFreq)
        linearBinCount = maxFreq / fftBinSize + 1;
    end
else
    linearBinCount = 0;
end

logBand_1 = bandCount;

% count the log space frequencies
lsFreq = centerFreq * lowSideMultiplier;
while (maxFreq > lsFreq)
    bandCount = bandCount + 1;
    logBinCount = logBinCount + 1;
    centerFreq =  getCenterFreq(base, bandsPerDivision, bandCount,firstOutputBandCenterFrequency);
    lsFreq  =  centerFreq * lowSideMultiplier;
end

freqTable = zeros((linearBinCount + logBinCount), 3);

% generate the linear frequencies
for i = 1:linearBinCount
    freqTable(i, 2)  =  bin1CenterFrequency + (i-1)*fftBinSize;
    freqTable(i, 1)  =  freqTable(i, 2) - fftBinSize/2;
    freqTable(i, 3)  =  freqTable(i, 2) + fftBinSize/2;
end

% generate the log spaced bands
for i = 1:logBinCount
    outBandNumber = linearBinCount + i;
    logBandNumber = logBand_1 + i - 1;
    freqTable(outBandNumber, 2) = getCenterFreq(base, bandsPerDivision,logBandNumber, firstOutputBandCenterFrequency);
    freqTable(outBandNumber, 1)  =  freqTable(outBandNumber, 2) * lowSideMultiplier;
    freqTable(outBandNumber, 3)  =  freqTable(outBandNumber, 2) * highSideMultiplier;
end

if (logBinCount > 0)
    % align the end of the linear with the start of the log:
    if (linearBinCount > 0)
        freqTable(linearBinCount, 3) = freqTable(linearBinCount+1, 1);
    end
    % handle end of data values:
    freqTable(outBandNumber, 3) = maxFreq;
    if (freqTable(outBandNumber, 2) > maxFreq)
        freqTable(outBandNumber, 2) = maxFreq;
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Utility function to compute the center frequencies with and without offset
% by half a band.
function centerFreq = getCenterFreq(base, bandsPerDivision, bandNum,firstOutBandCenterFreq)

if (bandsPerDivision == 10 || mod(bandsPerDivision, 2) == 1)
    centerFreq = firstOutBandCenterFreq * power(base, (bandNum-1) / bandsPerDivision);
else
    % from IEC 2014:
    G = power(10, .3);
    b = bandsPerDivision * 0.3;
    centerFreq = base*power(G, (2*(bandNum-1)+1)/(2*b));
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bandsOut = getBandSquaredSoundPressure(linLevel, fftBinSize, bin1CenterFrequency, firstBand, lastBand, freqTable )
% getBandSquaredSoundPressure - this function sums squared sound pressures
% to determine the in-band totals. The band edges are normally obtained
% from a call to 'getBandTable.m'
%
% inputs: linLevel - array of squared pressures from an FFT with a frequency
%		step size; index 1 (rows) are time, index 2 (columns) are freq.
%       fftBinSize - the size of the FFT bins in Hz.
%       bin1CenterFrequency: the freq in Hz of the first element of the FFT
%		array - normally this is frquency zero.
%       firstBand: the index in 'freqtable' of the first band to compute and
%		output
%       lastBand: the index in 'freqTable' of the last band to ocmpute and
%           output
%       freqTable - the list of band edges - Nx3 array where column 1 is the
%           lowest band frquency, column 2 is the center frequency and 3 is
%		the maximum.
%
% Outputs:  band squared sound pressure array with the same number of rows as
%	 linLevel and one column per band.
%
% Bruce Martin, JASCO Applied Sciences, Feb 2020.

nRows = size(linLevel, 1);

bandsOut = zeros(nRows, lastBand-firstBand+1);
step = fftBinSize / 2;
nFFTBins = size(linLevel, 2);
startOffset = floor(bin1CenterFrequency / fftBinSize);

for row = 1:nRows
    for j = firstBand:lastBand
        minFFTBin = floor((freqTable(j,1) / fftBinSize) + step) + 1 - startOffset;
        maxFFTBin = floor((freqTable(j,3) / fftBinSize) + step) + 1 - startOffset;
        if (maxFFTBin > nFFTBins)
            maxFFTBin = nFFTBins;
        end
        if (minFFTBin < 1)
            minFFTBin = 1;
        end

        if (minFFTBin == maxFFTBin)
            bandsOut(row, j) = linLevel(row, minFFTBin) *((freqTable(j, 3)- freqTable(j, 1))/ fftBinSize);
        else
            % Add the first partial FFT bin - take the top of the bin and
            % subtract the lower freq to get the amount we will use:
            % the top freq of a bin is bin# * step size - binSize/2 since bin
            %
            lowerFactor =((minFFTBin - step) * fftBinSize - freqTable(j, 1));
            bandsOut(row, j) = linLevel(row, minFFTBin) * lowerFactor;

            % Add the last partial FFT bin.
            upperFactor = freqTable(j, 3) - (maxFFTBin - 1.5*fftBinSize)* fftBinSize;
            bandsOut(row, j) = bandsOut(row, j) + linLevel(row, maxFFTBin)* upperFactor;
            %
            % Add any FFT bins in between min and max.
            if (maxFFTBin - minFFTBin) > 1
                bandsOut(row, j) = bandsOut(row, j) + sum( linLevel(row,minFFTBin+1:maxFFTBin-1) );
            end
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bandsOut = getBandMeanPowerSpectralDensity(linLevel, fftBinSize, bin1CenterFrequency, ...
    firstBand, lastBand, freqTable )
% getBandMeanPowerSpectralDensity - this function sums squared sound
%	pressures to determine the in-band totals then divides by the
%     bandwidths to get PSD. The band edges are normally obtained from a call
%	to 'getBandTable.m'
% getBandSquaredSoundPressure is called to get the band SPLs
%
% Note that the output of getBandSquaredSoundPressurere should satisfy
% Parseval's theorem, but the output of getBandMeanPowerSpectralDensity
% will not unless the bands are re-multiplied by the bandwidths.
% results are returned as linear units not levels.
%
% inputs: linLevel - array of squared pressures from an FFT with a frequency
%		step size; index 1 (rows) are time, index 2 (columns) are freq.
%       fftBinSize - the size of the FFT bins in Hz.
%       bin1CenterFrequency: the freq in Hz of the first element of the FFT
%		array  normally this is frquency zero.
%       firstBand: the index in 'freqtable' of the first band to compute and
%		output
%       lastBand: the index in 'freqTable' of the last band to ocmpute and
%           output
%       freqTable - the list of band edges - Nx3 array where column 1 is the
%           lowest band frquency, column 2 is the center frequency and 3 is
%		the maximum.
%
% Outputs:  band mean PSD array with the same number of rows as linLevel and
%       one column per band.

bandsOut = getBandSquaredSoundPressure(linLevel, fftBinSize, bin1CenterFrequency, firstBand, lastBand, freqTable );
nRows = size(linLevel, 1);
bandWidths = freqTable(:, 3) - freqTable(:, 1);
for row = 1:nRows
    bandsOut(row, :) = bandsOut(row, :) ./ bandWidths';
end
end