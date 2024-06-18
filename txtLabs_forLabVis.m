labDir = 'E:\S1099Hokkaido01_S01_RH403_WAV\SPICE_detector\clusters_below80K_3\cc10\TPWS_labels';
saveDir = 'E:\S1099Hokkaido01_S01_RH403_WAV\SPICE_detector\clusters_below80K_3\cc10';
saveName = 'S1099_below80K_Odontoceti_LOSpecific_labels';

labList = dir(fullfile(labDir,'*.mat'));

labelTable = [];

for i = 1:length(labList)

    load(fullfile(labDir,labList(i).name));
    zID = sortrows(zID,1);
%     labels = cellfun(@(x) erase(x,'Cluster'),labels,'UniformOutput',false);
    newAnnot = cell(size(zID,1),1);
    newAnnot(:,1) = {0};

    for k = 1:length(labels) % Get labels corresponding to cluster numbers
        lab = labels{1,k};
        idx = find(zID(:,2)==k);
        newAnnot(idx) = cellstr(lab);
    end

    clickEnd = zID(:,1) + datenum(0,0,0,0,0,0.001);
    newDat = table(zID(:,1),clickEnd,newAnnot);
    labelTable = [labelTable;newDat];

end

writetable(labelTable,fullfile(saveDir,[saveName,'.txt']));
