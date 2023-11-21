%% Calculate proportion of clicks included in cluster_bins mean spectra,
% and proportion isolated, across entire deployment

%directory containing cluster_bins output
inDir = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\ShortStur\Clusters\SS02\ClusterBins_w_NoiseDets\It_2';
binList = dir(fullfile(inDir,'*clusters*.mat'));

percSpec_perFile = zeros(length(binList),1,1);

for iA = 1:length(binList)
    
     load(fullfile(binList(iA).folder,binList(iA).name));
     percSpec = sum(horzcat(binData.nSpec))/sum(horzcat(binData.cInt));
     percSpec_perFile(iA) = percSpec;
    
end

figure
histogram(percSpec_perFile,0:0.01:1);
title('Proportion of Clicks Included in Clusters, Per File');
xlabel('%');
ylabel('Counts');

total_percSpec = mean(percSpec_perFile)*100;
total_percIso = (100-total_percSpec);

fprintf('%.2f %% of clicks included in mean spectra, %.2f %% isolated\n',...
    total_percSpec,total_percIso);