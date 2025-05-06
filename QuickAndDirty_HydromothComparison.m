
selDir = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\ChestnutSAVSoundscapes\HydromothComparison\SelectionTables';
selList = dir(fullfile(selDir,'*.txt'));
tones = [50,100,200,500,1000,2000,5000,10000,20000]';
bitDepth = 16;

compTable = [];
for i=1:length(selList)
    tb = readtable(fullfile(selDir,selList(i).name),'Delimiter','\t','VariableNamingRule','preserve');
    wavs = find(strcmp(table2array(tb(:,'View')),'Waveform 1'));
    tb = tb(wavs,:);

    reps = size(tb,1)/length(tones);
    tb.Freq = repmat(tones,reps,1);
    
    ind = strfind(selList(i).name,'_');
    ind = ind(end);
    tb.Name = cellstr(repmat(extractBefore(selList(i).name,ind),size(tb,1),1));

    if i==1
        bigTab = tb(:,{'Name','Freq','RMS Amp (U)'});
    else
        bigTab = [bigTab;tb(:,{'Name','Freq','RMS Amp (U)'})];
    end

end

nameMat = vertcat(bigTab.Name);

% Plot inter-device variability in sensitivity across freqs
figure(111),clf
boxplot(bigTab.("RMS Amp (U)")./((2^bitDepth)/2),bigTab.Name)
xticklabels({'Chestnut','ChestnutEdge','SAV1','SAV2'});
ylabel('Amplitude (% Full Scale)');
set(gca,'FontSize',14)


% For each freq, plot inter-device sensitivity
n=3;
m = ceil(length(tones)/n);
figure(222),clf
for i=1:length(tones)

    ind = find(bigTab.Freq==tones(i));

    subplot(n,m,i)
    boxplot(table2array(bigTab(ind,"RMS Amp (U)"))./((2^bitDepth)/2),bigTab.Name(ind))
    xticklabels({'Chestnut','ChestnutEdge','SAV1','SAV2'});
    ylabel('Amplitude (% Full Scale)');
    title(['Freq: ',num2str(tones(i)),' Hz'])
    set(gca,'FontSize',10)

end
