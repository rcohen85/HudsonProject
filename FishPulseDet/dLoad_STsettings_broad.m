function parametersST = dLoad_STsettings_broad
% Assign short term detector settings

parametersST.buff = 0.0025; % # of buffer in seconds to add on either side of area of interest
parametersST.chan = 1; % which channel do you want to look at?

parametersST.bpRanges = [5000 100000]; % Bandpass filter parameters in Hz [min,max]
parametersST.filterOrder = 5; % butterworth filter order used for band pass
parametersST.dBpp = 118; % minimum amplitude threshold in dB. 
parametersST.countThresh = []; % For predictability, keep this consistent between low and hi res steps.

parametersST.frameLengthUs = 2000; % For fft computation

%parametersST.frameLengthSec = .01; %Used for calculating fft size
parametersST.overlap = .50; % fft overlap
parametersST.REWavExt = '(\.x)?\.wav';%  expression to match .wav or .x.wav

% if you're using wav files that have a time stamp in the name, put a
% regular expression for extracting that here:
parametersST.DateRE = '_(\d*)_(\d*)';

% Examples:
%      a) File name looks like: "myFile_20110901_234905.wav" 
%         ie.:  "*_yyyymmdd_HHMMSS.wav"
%         So use:
%         parametersST.DateRE = '_(\d*)_(\d*)';
%
%      b) File name looks like: "palmyra102006-061104-012711_4.wav" 
%         ie.:  "*yyyy-yymmdd-HHMMSS*.wav"
%         So use:
%         parametersST.DateRE = '(\d{4})-\d{2}(\d{4})-(\d{6})';
% 

%%%%% GUIDED DETECTIONS? %%%%
parametersST.guidedDetector = 0; % flag to 1 if guided

