# Density Estimation Detector

Set of scripts to detect odontocete clicks above a threshold recieved level in .x.wav data. This two-pass detector is based on and incorporates previous work by Marie A. Roch, Simone Baumann-Pickering, and Sean M. Wiggins. 

Detector by Kait E. Frasier.

https://zenodo.org/badge/latestdoi/22897395

## How To Use This Code:

1. Place /DE_Detector/ and all subdirectories in Matlab path

2. Edit low resolution detector parameters 
		dLoad_STSettings.m

3. Edit high resolution detector parameters 
		dLoad_HRSettings.m

4. Edit the following lines at the top of *de_detector.m* to reflect your file locations, and which routines to run*

```matlab
% Set transfer function location
tfFullFile = 'E:\Code\TF_files\Recalculated\tf files\740_140303\740_140303_invSensit.tf';
% Note, if you don't have a transfer function just use:
% tfFullFile = [];

% Location of base directory containing directories of files to be analyzed
baseDir = 'H:\';

% Optional output directory location. Metadata directory will be created in outDir
% if specified, otherwise it will be created in baseDir.
% outDir = '<your path here>';
outDir  = 'I:\DCL\WAT_NC_'; 

% Name of the deployment. This should be the first few characters in the 
% directory(ies) you want to look in you want to look at. For now, 
% directory hierarchy is expected to be: basedir>depl*>*.x.wav
depl = 'WAT_NC_';

% Set flags indicating which routines to run. 
lowResDet = 1; %run short time detector.
highResDet = 1; %run high res detector

```

* NOTE: The high res detector relies on output from the low res step, but once you've run the low res, you don't need to re-run it, unless you change your parameters.

5. If you want to restrict the detector to only look at certain time periods, you can do so using a spreadsheet of start and end times specified in de_detector.m 

```matlab
%%%% Optional: guided detection spreadsheet, can be empty
gDxls = 'E:\Data\John Reports\DCLDEdata\WAT_NC_guidedDets.xlsx';
% gDxls = []; % if not used
```

## For .wav files

If you run the detector on directories of wav files, it will look for file start time information in the file name.
 
Edit the regular expression in the load settings scripts:

parametersHR.DateRE = '_(\d*)_(\d*)';

to match your filename date format. 

The result should be a string of numbers in the following order:
yyyymmddHHMMSS

The detector will determine what file type you're using, but to be on the safe side, I'd suggest analyzing wav and xwav files separately.


## Outputs

/metadata/<disk name> 
	contains all of the file types below:
- .c  - 
Low res detector output. This is a text file with flagged start and end times listed as 2 columns. Times are in seconds relative to .x.wav start time.

- .ctg  -   
High res detector output. This is a text file with detection start and end times listed as 2 columns. Times are in seconds relative to .x.wav start time.

- .ptg  - 
same as .ctg but has been run through a post processing step to remove redundant detections, etc.

- .mat  - 
matlab file containing various parameters describing the detected signals retained by all detection steps.

Note: Filenames match the name of the .xwav they describe.


