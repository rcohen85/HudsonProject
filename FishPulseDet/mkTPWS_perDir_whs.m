silbido_init
% mkTPWS_perDir.m

% Script takes output from de_detector.m and put it into a format for use in
% detEdit.
% Output:
% A TPWS.m file containing 4 variables:
%   MTT: An Nx2 vector of detection start and end times, where N is the
%   number of detections
%   MPP: An Nx1 vector of recieved level (RL) amplitudes.
%   MSP: An NxF vector of detection spectra, where F is dictated by the
%   parameters of the fft used to generate the spectra and any
%   normalization preferences.
%   f = An Fx1 frequency vector associated with MSP

clearvars

% Setup variables:
baseDir = 'I:\HAT_whs_over5\'; % directory containing de_detector output
outDir = 'I:\HAT_whs_over5\TPWS'; % directory where you want to save your TTPP file
% siteName = 'HAT03A'; % site name, only used to name the output file
% ppThresh = 120; % minimum RL in dBpp. If detections have RL below this
% threshold, they will be excluded from the output file. Useful if you have
% an unmanageable number of detections.

f = 5:.5:25; % freq band that the whistles were detected in, kHz

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check if the output file exists, if not, make it
if ~exist(outDir,'dir')
    fprintf('Creating output directory %s',outDir)
    mkdir(outDir)
end
%letterCode = 97:122;
dirSet = dir(fullfile(baseDir,'*_disk*'));
for itr0 = 1:length(dirSet)
    if dirSet(itr0).isdir &&~strcmp(dirSet(itr0).name,'.')&&...
            ~strcmp(dirSet(itr0).name,'..')
        
        %letterFlag = 0; % flag for knowing if a letter should be appended to disk name
        inDir = fullfile(baseDir,dirSet(itr0).name);
        fileSet = what(inDir);
        lfs = length(fileSet.mat);
        timesVec = [];
        snrVec = [];
        specVec = [];
        contourVec = [];
        for itr2 = 1:lfs
            thisFile = fileSet.mat(itr2);
            
            load(char(fullfile(inDir,thisFile)),'-mat','tonals','PARAMS')
            nTonals = tonals.size();
            snrTonal = zeros(nTonals,1);
            timeTonal = zeros(nTonals,1);
            specTonal = zeros(nTonals, length(f));
            contourTonal = nan(nTonals, 200);

            for iT = 0:nTonals-1
                thisTonal = tonals.get(iT);
                tonalFreq = thisTonal.get_freq();
                snrTonal(iT+1,1) = mean(thisTonal.get_snr());
                timeTonal(iT+1,1) = mean(thisTonal.get_time());
                specTonal(iT+1,:) = histc(tonalFreq/1000,f);
                tTrunc = min(200,length(tonalFreq));
                contourTonal(iT+1,1:tTrunc) = tonalFreq(1:tTrunc);
            end
            
            if nTonals>0
                % fileStart = datenum(hdr.start.dvec);
                % HAT03A_DL34_130611_190615
                tfc = char(thisFile);
                fileStart = datenum([2000+str2num(tfc(13:14)),str2num(tfc(15:16)),...
                    str2num(tfc(17:18)),str2num(tfc(20:21)),str2num(tfc(22:23)),...
                    str2num(tfc(24:25))]);
                posDnum = (timeTonal/(60*60*24)) + fileStart;
                
                timesVec = [timesVec; posDnum];
                snrVec = [snrVec; snrTonal];
                specVec = [specVec; specTonal];
                contourVec = [contourVec;contourTonal];
                
            end
            fprintf('Done with file %d of %d \n',itr2,lfs)
        end
        
        MSN = contourVec;
        MTT = timesVec;
        MPP = snrVec;
        MSP = specVec;
        ttppOutName =  [fullfile(outDir,dirSet(itr0).name),'_whs_TPWS1','.mat'];
        
        save(ttppOutName,'MTT','MPP','MSP','MSN','f','-v7.3')
        
        MTT = [];
        MPP = [];
        MSP = [];
        MSN = [];
        
        timesVec = [];
        snrVec = [];
        specVec = [];
        contourVec = [];
        
        fprintf('Done with directory %d of %d \n',itr0,length(dirSet))
        
    end
end


