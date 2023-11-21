
TPWSdir = 'E:\S1099Hokkaido01_S01_RH403_WAV\SPICE_detector\TPWSbelow80K';
labDir = 'E:\S1099Hokkaido01_S01_RH403_WAV\SPICE_detector\clusters_below80K\cc';

TPWSlist = dir(fullfile(TPWSdir,'*TPWS1.mat'));
labList = dir(fullfile(labDir,"*ID1.mat"));

allPulses = struct("MSP",[],"MTT",[],"MPP",[],"MSN",[]);
stats_perFile = zeros(length(TPWSlist),3,1);

for i = 1:length(TPWSlist)

    load(fullfile(TPWSdir,TPWSlist(i).name));
    load(fullfile(labDir,labList(i).name));
    zID = sortrows(zID,1);
    labels = cellfun(@(x) erase(x,'Cluster'),labels,'UniformOutput',false);

    stats_perFile(i,1) = size(zID,1);
    stats_perFile(i,2) = size(MTT,1);
    stats_perFile(i,3) = size(zID,1)/size(MTT,1);


    for j = 1:length(labels)
        lab = str2double(labels{1,j});
        idx = find(zID(:,2)==j);
        [Lia,Locb] = ismember(zID(idx,1),MTT);

        if ~isempty(idx)
            specs = MSP(Locb,:);
            specs = specs - min(specs,[],2);
            specs_norm = specs./max(specs,[],2);

            if lab > length(allPulses)
                allPulses(lab).MSP = []; % initialize rows as needed
            end
            
            allPulses(lab).MSP = [allPulses(lab).MSP;specs_norm];
            allPulses(lab).MTT = [allPulses(lab).MTT;MTT(Locb)];
            allPulses(lab).MPP = [allPulses(lab).MPP;MPP(Locb)];
            allPulses(lab).MSN = [allPulses(lab).MSN;MSN(Locb,:)];
            
        end

    end

    MSP = [];
    MSN = [];
    MPP = [];
    MTT = [];
    zID = [];

end

save(fullfile(labDir,'CompiledPulses.mat'),'allPulses','-v7.3');

% meanSpecs = [];
% for i = 1:length(labels)
%     meanSpecs(i,:) = mean(allPulses(i).MSP,1);
% end

% Plot concatenated pulses for each type
n = 4; % number of columns of subplots, one subplot per cluster
m = ceil(length(labels)/n); % number of rows of subplots
figure(999)
set(gcf,'Position',[300 50 1600 1200])
for i = 1:length(labels)
subplot(m,n,i)
[~,I] = max(allPulses(i).MSP,[],2);
[I,sortInd] = sort(I);
imagesc([],f,allPulses(i).MSP(sortInd,:)');set(gca,'YDir','Normal')
xlabel('Pulse Number')
ylabel('Frequency (kHz)');
title(sprintf('Cluster %d',i));
end
saveas(gcf,fullfile(labDir,'CompiledPulses.png'));

% Plot percent of labeled pulses per file
total_percLabeled = (sum(stats_perFile(:,1))/sum(stats_perFile(:,2)))*100;
total_percIso = (100-total_percLabeled);

figure(99)
histogram(stats_perFile(:,3),0:0.01:1);
title('Percent of Pulses with Labels, Per File');
xlabel('%');
ylabel('Counts');
text(0,2,[num2str(round(total_percLabeled,2)),'% of pulses labeled, ',num2str(round(total_percIso,2)),'% isolated']);
saveas(gcf,fullfile(labDir,'PercentLabeled.png'));
