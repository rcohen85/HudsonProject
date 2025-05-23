function parametersST = dLoad_STsettings_FishPulse
% Assign short term detector settings

parametersST.buff = 0.1; % # of buffer in seconds to add on either side of area of interest
parametersST.chan = 1; % which channel do you want to look at?

parametersST.bpRanges = [10 10000]; % Bandpass filter parameters in Hz [min,max]
parametersST.filterOrder = 5; % butterworth filter order used for band pass
parametersST.dBpp = []; % minimum amplitude threshold in dB. 
parametersST.countThresh = 40; % For predictability, keep this consistent between low and hi res steps.

parametersST.frameLengthUs = 50000; % For fft computation
parametersST.mergeThr = 1000;% min gap between energy peaks in us

parametersST.frameLengthSec = .05; %Used for calculating fft size
parametersST.overlap = .50; % fft overlap
parametersST.REFileExt = '.wav';%  expression to match file extension

% if you're using files that have a time stamp in the name, put a
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

