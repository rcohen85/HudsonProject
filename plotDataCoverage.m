% inDirs = {'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\AtlSturSpw\ST600\AS09';
%     'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\AtlSturSpw\ST600\AS10';
%     'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\AtlSturSpw\ST600\AS11';
%     'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\AtlSturSpw\ST600\AS12';
%     'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\BlackCreek\Swift\BC10';
%     'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\BlackCreek\Swift\BC11';
%     'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\BlackCreek\Swift\BC12';
%     'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\BlackCreek\Swift\BC13';
%     'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\Iona\ST600\ICh01';
%     'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\Iona\ST600\ICh02';
%     'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\Iona\ST600\ICh03';
%     'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\Iona\ST600\ICh04';
%     'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\Iona\Swift\II01';
%     'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\Iona\Swift\II02';
%     'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\Iona\Swift\II03';
%     'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\Piermont\ST600\PCh01';
%     'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\Piermont\ST600\PCh02';
%     'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\Piermont\ST600\PCh03';
%     'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\Piermont\ST600\PCh04';
%     'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\Piermont\ST600\PCh05';
%     'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\Piermont\Swift\PM01';
%     'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\Piermont\Swift\PM02';
%     'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\Piermont\Swift\PM03';
%     'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\Stockport\ST600\SCh01';
%     'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\Stockport\ST600\SCh02';
%     'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\Stockport\Swift\SCr01';
%     'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\Stockport\Swift\SCr02';
%     'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\Stockport\Swift\SCr03';
%     'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\Tivoli\ST600\TCh01';
%     'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\Tivoli\ST600\TCh02';
%     'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\Tivoli\ST600\TCh04';
%     'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\Tivoli\Swift\TBN02';
%     'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\Tivoli\Swift\TBN03';
%     'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\Tivoli\Swift\TBS02';
%     'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\Tivoli\Swift\TBS03';
%     'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\Tivoli\Swift\TBN04';
%     'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\Tivoli\Swift\TBS04'};
% inDirs = {'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_FLAC\AtlSturSpw\ST300\AS02';
%     'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_FLAC\AtlSturSpw\ST300\AS03';
%     'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_FLAC\AtlSturSpw\ST300\AS04';
%     'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_FLAC\AtlSturSpw\ST300\AS05'};
% inDirs = {'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\NorriePoint\Hydromoth\NP13';
%     'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\NorriePoint\Hydromoth\NP12';
%     'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\NorriePoint\Hydromoth\NP11';
%     'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\NorriePoint\Hydromoth\NP10';
%     'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\NorriePoint\Hydromoth\NP09'};

% inDirs = {'W:\projects\2024_NPS-GLBA_AK_173170\Sounds\Terrestrial\2024 GLBALSTL Lower South Tidal Inlet\01 DATA\AUDIO';
%     'W:\projects\2024_NPS-GLBA_AK_173170\Sounds\Terrestrial\2024 GLBAUSTL
%     Upper South Tidal Inlet\01 DATA\AUDIO'};
inDirs = {'W:\projects\2024_NPS-GLBA_AK_173170\Sounds\Terrestrial\2024 GLBATARR Tarr Inlet\01 DATA\AUDIO'};

dates = [];
durs = [];
names = [];
for i=1:size(inDirs,1)
    
    ind = strfind(inDirs{i},'\');
    ind = ind(end);
    siteName = extractAfter(inDirs{i},ind);
%     siteName = extractBefore(siteName,length(siteName)-1);
    fileList = dir(fullfile(inDirs{i},'*.MP3'));
    fileNames = cellstr(vertcat(fileList.name));
    for j=1:size(fileList,1)
        info = audioinfo(fullfile(fileList(j).folder,fileList(j).name));
        durs = [durs;vertcat(info.TotalSamples)./vertcat(info.SampleRate)];
    end
    fileStarts = cellfun(@(x) regexp(x,'_\d{8}[_]\d{6}','match'),fileNames,'UniformOutput',0);
    fileStarts = (datetime([fileStarts{:}],'InputFormat','_yyyyMMdd_HHmmss'))';
    dates = [dates;fileStarts];
    names = [names;cellstr(repmat(siteName,size(fileStarts,1),1))];
    fprintf('Done with directory %d of %d\n',i,size(inDirs,1))

end

% nameCat = findgroups(names);
nameCat = NaN(length(names),1);
nameCat(cell2mat(strfind(names,'SCh'))) = 11;

% TO DO: Group by site name before summing duration
totalDur = sum(durs)/(60*60*24);

figure(3),clf
% plot(dates,repmat(i,length(dates),1),'.','MarkerSize',20)
plot(dates,nameCat,'.','MarkerSize',20)
xline([datetime(2023,07,08),datetime(2023,07,25)],'LineWidth',2,'Color','r')
yticks([unique(nameCat)]);
% yticklabels(unique(names));
yticklabels({'WC1','WC2','Open Water','SAV1','SAV2'});
ylim([min(nameCat)-1,max(nameCat)+1])
% xlim([datetime(2023,03,01),datetime(2024,08,15)])
% xticks([datetime(2021,04,15),datetime(2021,05,01),datetime(2021,05,15),datetime(2021,06,01), ...
%     datetime(2021,06,15),datetime(2021,07,01),datetime(2021,07,15),datetime(2021,07,30)])
set(gca,'FontSize',18)
