function fishPulse_detector
% This is the starting point for a simplified detector based on Marie Roch's  
% teager energy detector. It includes ideas/code snips from Simone, 
% and calls functions from triton. 
% The goal of this detector is to have predictable performance, 
% for use with model-based density estimation efforts. To accomplish this,
% it uses a simple energy threshold to identify clicks, thereby reducing
% the impact of changing noise conditions on detectability. 

% Known issue: Prop cavitation noise often makes it through detector and
% classifier steps.

% The low and hi-res detection passes still happen, but no teager energy 
% is used.

% All input parameters are contained within two separate scripts:
%   dLoad_STsettings : settings for low res detector
%   dLoad_HRsettings : settings for hi res detector
% See those files for info on settings.

clearvars
fclose all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set transfer function location
% tfFullFile = 'G:\Shared drives\MBARC_TF\600-699\673\673_120705_A_HARP.tf';
% Note, if you don't have a transfer function just use:
tfFullFile = '';%176.3;
deviceType = [];%'ST';

% Location of base directory containing directories of files to be analyzed
baseDir = 'C:\Users\rec297\Desktop\Test\Wav';

% Optional output directory location. Metadata directory will be created in outDir
% if specified, otherwise it will be created in baseDir.
% outDir = '<your path here>';
outDir  = 'C:\Users\rec297\Desktop\Test\Wav\S1099'; 

[metaDir,~] = dBuild_dirs(baseDir,outDir);
% inDisk = fileparts(baseDir(1:3));
% diary(fullfile(metaDir,sprintf('diary_%s.txt',datestr(now,'YYYYMMDD'))))

% Name of the deployment. This should be the first few characters in the 
% directory(ies) you want to look in you want to look at. For now, 
% directory hierarchy is expected to be: basedir>depl*>*.x.wav
% TODO: implement recursive directory search for more flexibility.
depl = 'S1099';

% Set flags indicating which routines to run. 
lowResDet = 1; %run short time detector.
highResDet = 1; %run high res detector

%%%% Optional: guided detection spreadsheet, can be empty
% gDxls = 'E:\Data\John Reports\DCLDEdata\WAT_NC_guidedDets.xlsx';
gDxls = []; % if not used

%%%%%%%%%%%%%%%%%%%%%%%%%%%% End Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
paramsST = [];
paramsHR = [];
gD = 0;
% load settings
if lowResDet
    paramsST = dLoad_STsettings_FishPulse;
    if paramsST.guidedDetector
       gD = 1; % if guided detector is true in either place, 
       % and the stage (low or high res) it's in is activated, then it
       % needs to be true any other active stage.
    end
end
if highResDet
    paramsHR = dLoad_HRsettings_FishPulse;
    if paramsHR.guidedDetector
       gD = 1;
    end
end


% Build list of (x)wav names in the base directory.
% Right now only wav and xwav files are looked for.
if lowResDet
    detFiles = dFind_files(baseDir,depl,metaDir,paramsST.REFileExt);
elseif highResDet
    detFiles = dFind_files(baseDir,depl,metaDir,paramsHR.REFileExt);
end

if gD && ~isempty(gDxls)
    [detFiles,encounterTimes] = guidedDetection(detFiles,gDxls);
    fprintf('Using guided detections from file %s \n',gDxls')
    %graphDir = 1;
else 
    encounterTimes = [];
end

% return a list of files to be built
if lowResDet
    [fullFiles,fullLabels] = get_fileset(baseDir,metaDir,detFiles,paramsST.REFileExt);
elseif highResDet
    [fullFiles,fullLabels] = get_fileset(baseDir,metaDir,detFiles,paramsHR.REFileExt); 
end

% profile on
% profile clear
if ~isempty(detFiles)
    % Short time detector
    if lowResDet
        tic 
        display('Beginning low-res detection')
        dtST_batch(fullLabels,fullFiles,paramsST,tfFullFile,deviceType); % run detector
        display('Done with low-res detector')
        toc
    end
    
    % High res detector
    if highResDet
        tic
        display('Beginning high-res detection')
        dHighres_click_batch(fullFiles,fullLabels,paramsHR,...
            tfFullFile,encounterTimes)
        display('Done with high-res detector')
        toc
    end
else
    disp('Error: No wav/xwav files found')
end

% profile viewer
% profile off
% diary('off')