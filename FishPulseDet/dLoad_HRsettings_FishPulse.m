function parametersHR = dLoad_HRsettings_broad

%%% BANDPASS FILTER AND FFT PARAMS %%
parametersHR.bpRanges = [10,10000]; % Bandpass filter params in Hz [min,max]
parametersHR.filterOrder = 5; % butterworth filter order used for band pass
parametersHR.frameLengthUs = 50000; % For fft computation
parametersHR.overlap = .5; % FFT overlap (in decimal, not percent form)
parametersHR.chan = 1; % which channel do you want to look at?
parametersHR.clipThreshold = .98;%  Normalized clipping threshold btwn 0 and 1.  If empty, 
% assumes no clipping. 

%%% Saving options %%%
parametersHR.saveNoise = 0; % Make 1 if you want to save noise samples with each click. 
% Beware: this can make big files if you have a lot of detections.
parametersHR.saveForTPWS = 1; % Save just enough data to build TPWS files. Should help
% limit metadata size.

%%% RECIEVED LEVEL THRESHOLD %%%
parametersHR.dBpp = 110;% minimum  RL threshold - dB peak to peak.
parametersHR.countThresh = [];%40; 
%%% CLICK ENVELOPE PARAMS %%%
parametersHR.energyThr = 0.7; % n-percent energy threshold for high-energy envelope duration
parametersHR.dEvLims = [];  % [min,max] Envelope energy distribution comparing 
% first half to second half of high energy envelope of click. If there is
% more energy in the first half of the click (dolphin) dEv >0, If it's more
% in the second half (boats?) dEv<0. If it's about the same (beaked whale)
% dEnv ~= 0 , but still allow a range...
parametersHR.envDurLims = [100,120000];% [min,max] duration in microsec 
% allowed for high energy envelope of click


%%% OTHER DETECTION THRESHOLD PARAMS %%%
% parametersHR.bw3dbMin = 2; % In kHz. discard click if less than this wide. 
% Helps with false detections from things like echosounders and rain.
parametersHR.cutPeakBelowKHz = 0.1; % discard click if peak frequency below X kHz
parametersHR.cutPeakAboveKHz = 10; % discard click if peak frequency above Y kHz 
parametersHR.minClick_us = 5000;% Minimum duration of a click in us 
parametersHR.maxClick_us = 160000; % Max duration of a click including echos
parametersHR.maxNeighbor = 3600; % max time in seconds allowed between neighboring clicks

parametersHR.mergeThr = 1000;% min gap between energy peaks in us. Anything less
% will be merged into one detection the beginning of the next is fewer
% samples than this, the signals will be merged.

parametersHR.energyPrctile = 70; % sets the threshold at which click start 
% and end are defined. In a time window of interest, the detector finds the
% maximum energy peak, and then walks forward and back in time along a smoothed
% energy vector, until it gets to an energy level set by this threshold, 
% ie. 70th percentile of energy in the time window. 

%%% GUIDED DETECTOR OPTIONS %%%%%%%
parametersHR.guidedDetector = 0; % set to 1 for guided detections.

%%% FILE EXTENSION %%%%%%
% need this in here in case short time detector is not run
parametersHR.REFileExt = '.flac';%  expression to match file extension

%%% FILE NAMING OPTION %%%%%%%%
% if you're using files that have a time stamp in the name, put a
% regular expression for extracting that here:
parametersHR.DateRE = '_(\d*)_(\d*)';

% mine look like "filename_20110901_234905.wav" 
% ie "*_yyyymmdd_HHMMSS.wav"


%%% POST PROCESSING FLAGS %%%%%%%%
% 1 for true, 0 for false
parametersHR.rmLonerClicks = 0;
parametersHR.rmEchos = 0;
if parametersHR.rmEchos
    parametersHR.lockOut = 0.05; %min gap between clicks
end
%%% Output file extensions. Probably don't need to be changed %%%
parametersHR.clickAnnotExt = 'cTg';
parametersHR.ppExt = 'pTg';
parametersHR.groupAnnotExt = 'gTg';
