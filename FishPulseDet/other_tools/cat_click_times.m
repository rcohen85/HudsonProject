% cat_click_times.m
% Could turn this into a function and instert at the end of de_detector.m
% function cat_click_times(inDir)
inDir = 'E:\metadata\bigDL'; % the path to your directory of detector outputs goes here
matList = dir(fullfile(inDir,'Tin*.mat')); % Add wildcard to match the files you want to process.

clickDnum = [];
sec2dnum = 60*60*24; % conversion factor to get from seconds to matlab datenum

% iterate over detector-derived mat files in directory
for i1 = 1:length(matList)
    clickTimes = [];
    clickDnumTemp = [];
    % only need to load hdr and click times
    load(fullfile(inDir,matList(i1).name),'hdr','clickTimes')
    
    
    if ~isempty(clickTimes)
        % determine true click times 
        clickDnumTemp = (clickTimes./sec2dnum) + hdr.start.dnum + datenum([2000,0,0]);
        clickDnum = [clickDnum,clickDnumTemp]; %save to one vector
        
        % write label file:
        clickTimeRel = zeros(size(clickDnumTemp));
        rawStarts = hdr.raw.dnumStart + datenum([2000,0,0]);
        rawDurs = (hdr.raw.dnumEnd-hdr.raw.dnumStart)*sec2dnum;
        % generate label file by replacing .mat extension with .lab for
        % wavesurfer:
        outFileName = strrep(matList(i1).name,'.mat','.lab'); 
        
        % open file for writing
        fidOut = fopen(fullfile(inDir,outFileName),'w+');
        for i2 = 1:size(clickTimes,1)
            % figure out the closest raw start less than the click start,
            % and subtract that time out of the click time, so it's not
            % relative
            thisRaw = find(rawStarts<=clickDnumTemp(i2,1),1,'last');
            clickTimeRel = (clickDnumTemp(i2,:) - rawStarts(thisRaw))*sec2dnum;
            % add back in the duration of recordings in seconds preceeding
            % this raw file
            if thisRaw>1
                clickTimeRel = clickTimeRel + sum(rawDurs(1:thisRaw-1));
                
            end
            % writes click number as third item in row, because some label
            % is needed:
            fprintf(fidOut, '%f %f %d\n', clickTimeRel(1,1),clickTimeRel(1,2),i2);
            
        end
        fclose(fidOut);
    end
    
end

save(fullfile(inDir,'AllClickDnum.mat'),'clickDnum')