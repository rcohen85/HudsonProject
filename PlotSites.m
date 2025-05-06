%% Plot NERRS component sites as polygons on map
polyLats = {[42.36,42.36,42.28,42.28,42.36],... % Stockport
    [ 42.07, 42.07,41.99,41.99,42.07],... % Tivoli
    [ 41.88, 41.88,41.79,41.79, 41.88],... % Norrie Point
    [ 41.35,41.35,41.27,41.27,41.35],... % Iona
    [ 41.26,41.26,41.21, 41.21,41.26],... % Furnacebrook & Enderkill Creeks
    [41.07,41.07,40.99,40.99,41.07]}; % Piermont
polyLons = {[-73.85,-73.73,-73.73,-73.85,-73.85],... % Stockport
    [-73.99,-73.87,-73.87,-73.99,-73.99],... % Tivoli
    [-74.01,-73.89,-73.89,-74.01,-74.01],... % Norrie Point
    [-74.02,-73.90,-73.90,-74.02,-74.02],... % Iona
    [-73.96,-73.87,-73.87,-73.96,-73.96],... % Furnacebrook & Enderkill Creeks
    [-73.96,-73.84,-73.84,-73.96,-73.96]}; % Piermont
sitePoly = geopolyshape(polyLats,polyLons);

figure(66), clf
h = geoplot(sitePoly,'FaceColor','none','EdgeColor','b','LineWidth',2)
geobasemap topographic
set(gca,'FontSize',14)
t=h.Parent;
t.LatitudeLabel.String="";
t.LongitudeLabel.String="";
t.LongitudeAxis.TickLabelFormat = '-dd';
t.LatitudeAxis.TickLabelFormat = '-dd';

savename = 'P:\users\cohen_rebecca_rec297\CCB\HudsonProject\AllHudson_SitePolygons_topoMap_FBEK.png';
saveas(gcf,savename,'png');

%% Plot all NERRS Hudson River Valley sites (main channel & tributaries)
siteNames = {'Norrie Point','NP';'Atlantic Sturgeon Spawning','AS';'Black Creek','BC';'Enderkill Creek','EC';...
    'Furnace Brook, Above Dam','FBAD';'Furnace Brook, Below Dam','FBBD';'Shortnose Sturgeon Overwintering','SS';...
    'Tivoli Bay North','TBN';'Tivoli Bay South','TBS';'Tivoli Channel','TCh';'Iona Island','II';'Iona Channel','ICh';...
    'Piermont Marsh','PM';'Piermont Channel','PCh';'Stockport Creek','SCr';'Stockport Channel','SCh'};

deployInfo = readtable('W:\projects\2022_NOAA-NERRS_HudsonNY_144488\DeploymentInfo.xlsx','VariableNamingRule','preserve');
for i=1:height(deployInfo)
    if strcmp(deployInfo.Site(i),'')
        deployInfo.Site(i) = deployInfo.Site(i-1);
    end
end
deployInfo([27,53,79:80],:) = [];
% sites = unique(deployInfo.Site);
siteLats = [];
siteLons = [];
for i=1:length(siteNames)
%     siteInds = find(strcmp(sites{i},deployInfo.Site));
    siteInds = find(contains(deployInfo.Deployment,siteNames{i,2}));
    siteLats(i) = mean(deployInfo.Lat(siteInds),'omitnan');
    siteLons(i) = mean(deployInfo.Lon(siteInds),'omitnan');
end

figure(88),clf
h = geoscatter(siteLats,siteLons,60,'filled')
geobasemap topographic
t=h.Parent;
t.LatitudeLabel.String="";
t.LongitudeLabel.String="";

savename = 'P:\users\cohen_rebecca_rec297\CCB\HudsonProject\AllHudson_SiteMap_Larger.png';
saveas(gcf,savename,'png');

%% Plot Stockport sites
names = {'Stockport Channel','Stockport Creek'}';
devices = categorical({'ST600','Swift'}');
lats = [42.31825,42.309444]';
lons = [-73.775278,-73.771389]';
% siteTable = table(names,devices,lats,lons);
cmap = lines(length(unique(devices)));

figure(99),clf
h = geoscatter(lats,lons,200,cmap(devices,:),'filled')
geobasemap topographic
geolimits([42.29,42.34],[-73.815,-73.725])
% hold on % workaround to display discrete color legend instead of colorbar
% for i=1:size(cmap,1)
% h(i)=geoscatter(lats(1),lons(1),0.5,cmap(i,:),'filled');
% end
% legend(h, unique(devices))
% hold off
t=h.Parent;
t.LatitudeLabel.String="";
t.LatitudeAxis.FontSize=25;
t.LongitudeLabel.String="";
t.LongitudeAxis.FontSize=25;
t.LongitudeAxis.TickLabelFormat = '-dd';
t.LatitudeAxis.TickLabelFormat = '-dd';

savename = 'P:\users\cohen_rebecca_rec297\CCB\HudsonProject\StockportSites.png';
saveas(gcf,savename,'png');

%% Plot Tivoli sites
names = {'Tivoli Bay North','Tivoli Bay South','Tivoli Channel'}';
devices = categorical({'Swift','Swift','ST600'});
lats = [42.039596,42.0178886666667,42.020556]';
lons = [-73.9171333333333,-73.918852,-73.9295]';
% siteTable = table(names,devices,lats,lons);
cmap = lines(length(unique(devices)));

figure(99),clf
h = geoscatter(lats,lons,200,cmap(devices,:),'filled')
geobasemap topographic
geolimits([42,42.06],[-73.97,-73.88])
% hold on % workaround to display discrete color legend instead of colorbar
% for i=1:size(cmap,1)
% h(i)=geoscatter(lats(1),lons(1),0.5,cmap(i,:),'filled');
% end
% legend(h, unique(devices))
% hold off
t=h.Parent;
t.LatitudeLabel.String="";
t.LatitudeAxis.FontSize=25;
t.LongitudeLabel.String="";
t.LongitudeAxis.FontSize=25;
t.LongitudeAxis.TickLabelFormat = '-dd';
t.LatitudeAxis.TickLabelFormat = '-dd';

savename = 'P:\users\cohen_rebecca_rec297\CCB\HudsonProject\TivoliSites.png';
saveas(gcf,savename,'png');

%% Plot Norrie Point sites
% names = {'NorriePoint','ShortnoseStur','AtlanticStur','BlackCreek','Enderkill Creek'};
names = {'AtlanticStur','Telem'};
% devices = categorical({'Swift','ST600','ST600','Swift','Swift'});
devices = categorical({'ST300','TelemReceiver'});
% lats = [41.83167,41.86768,41.81440,41.82400,41.839]';
lats = [41.81440,41.83077063];
% lons = [-73.94194,-73.93382,-73.94540,-73.95900,-73.935]';
lons = [-73.94540,-73.94737395];
% siteTable = table(names,devices,lats,lons);
cmap = lines(length(unique(devices)));
% 
% figure(9)
% gb1 = geobubble(siteTable.lats,siteTable.lons,1);
% gb1.BubbleWidthRange = 1;
% gb1.Basemap = 'topographic';
% geolimits([41.8,41.8833],[-73.9750 -73.9083])

figure(9999),clf
h = geoscatter(lats,lons,100,cmap(devices,:),'filled')
geobasemap topographic
geolimits([41.79, 41.85],[-73.965,-73.935])
% hold on % workaround to display discrete color legend instead of colorbar
% for i=1:size(cmap,1)
% h(i)=geoscatter(lats(1),lons(1),0.5,cmap(i,:),'filled');
% end
% legend(h, unique(devices))
% hold off
t=h.Parent;
t.LatitudeLabel.String="";
t.LatitudeAxis.FontSize=14;
t.LongitudeLabel.String="";
t.LongitudeAxis.FontSize=14;
t.LongitudeAxis.TickLabelFormat = '-dd';
t.LatitudeAxis.TickLabelFormat = '-dd';

% gb2 = geobubble(siteTable.lats,siteTable.lons,1,siteTable.devices)
% gb2.BubbleWidthRange = 10;
% gb2.ColorLegendTitle = 'Device';
% gb2.Basemap = 'topographic';

savename = 'P:\users\cohen_rebecca_rec297\CCB\HudsonProject\NorriePointSites.png';
savename = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AtlStur\Figures\SiteMap.pdf'
saveas(gcf,savename,'png');

%% Plot Atlantic sturgeon site

% lats = [41.81440]';
% lons = [-73.94540]';
% polyLats = [ 41.818, 41.818,41.814086,41.814086,41.818];
% polyLons = [-73.9458,-73.9443,-73.9454,-73.9469,-73.9458];
% sitePoly = geopolyshape(polyLats,polyLons);
% 
% h = geoplot(sitePoly,'FaceColor','none','EdgeColor','w','LineWidth',2)
% hold on
% geoscatter(lats,lons,200,'filled')
% hold off
% geobasemap topographic
% geolimits([41.8098, 41.824],[-73.95,-73.939])
% set(gca,'FontSize',14)
% t=h.Parent;
% t.LatitudeLabel.String="";
% t.LongitudeLabel.String="";
% t.LongitudeAxis.TickLabelFormat = '-dd';
% t.LatitudeAxis.TickLabelFormat = '-dd';

polyLats = [ 41.8214, 41.8214,41.8124,41.8124,41.8214];
polyLons = [-73.9487,-73.9395,-73.9395,-73.9487,-73.9487];
sitePoly = geopolyshape(polyLats,polyLons);

figure(777),clf
h = geoplot(sitePoly,'FaceColor','none','EdgeColor','b','LineWidth',2)
geobasemap topographic
geolimits([41.8, 41.88],[-74,-73.89])
t=h.Parent;
t.LatitudeLabel.String="";
t.LatitudeAxis.FontSize=15;
t.LongitudeLabel.String="";
t.LongitudeAxis.FontSize=15;
t.LongitudeAxis.TickLabelFormat = '-dd';
t.LatitudeAxis.TickLabelFormat = '-dd';


savename = 'P:\users\cohen_rebecca_rec297\CCB\HudsonProject\AtlanticSturgeonSite';
saveas(gcf,savename,'png');
%% Plot Iona sites
names = {'Iona Island','Iona Channel'}';
devices = categorical({'Swift','ST600'});
lats = [41.303138,41.29806]';
lons = [-73.979325,-73.96639]';
% siteTable = table(names,devices,lats,lons);
cmap = lines(length(unique(devices)));

figure(99),clf
h = geoscatter(lats,lons,200,cmap(devices,:),'filled')
geobasemap topographic
geolimits([41.29,41.32],[-74.0,-73.94])
% hold on % workaround to display discrete color legend instead of colorbar
% for i=1:size(cmap,1)
% h(i)=geoscatter(lats(1),lons(1),0.5,cmap(i,:),'filled');
% end
% legend(h, unique(devices))
% hold off
t=h.Parent;
t.LatitudeLabel.String="";
t.LatitudeAxis.FontSize=33;
t.LongitudeLabel.String="";
t.LongitudeAxis.FontSize=33;
t.LatitudeAxis.TickValues = [41.3,41.31]
t.LongitudeAxis.TickLabelFormat = '-dd';
t.LatitudeAxis.TickLabelFormat = '-dd';

savename = 'P:\users\cohen_rebecca_rec297\CCB\HudsonProject\IonaSites2.png';
saveas(gcf,savename,'png');

%% Plot Furnacebrook Creek & Enderkill Creek
names = {'Enderkill Creek','Furnace Brook, Above Dam','Furnace Brook, Below Dam'};
devices = categorical({'Swift','Swift','Swift'});
lats = [41.839,41.231,41.22798]';
lons = [-73.935,-73.918,-73.920835]';
% siteTable = table(names,devices,lats,lons);
cmap = lines(2);

figure(99),clf
h = geoscatter(lats,lons,200,cmap(2,:),'filled')
geobasemap topographic
geolimits([41.22,41.24],[-73.96,-73.88])
% hold on % workaround to display discrete color legend instead of colorbar
% for i=1:size(cmap,1)
% h(i)=geoscatter(lats(1),lons(1),0.5,cmap(i,:),'filled');
% end
% legend(h, unique(devices))
% hold off
t=h.Parent;
t.LatitudeLabel.String="";
t.LatitudeAxis.FontSize=25;
t.LongitudeLabel.String="";
t.LongitudeAxis.FontSize=25;
t.LatitudeAxis.TickValues = [41.22,41.23,41.24]
t.LongitudeAxis.TickLabelFormat = '-dd';
t.LatitudeAxis.TickLabelFormat = '-dd';

savename = 'P:\users\cohen_rebecca_rec297\CCB\HudsonProject\FBEKSites.png';
saveas(gcf,savename,'png');

%% Plot Piermont Sites
names = {'Piermont Marsh','Piermont Channel'}';
devices = categorical({'Swift','ST600'});
lats = [41.030574,41.0357]';
lons = [-73.911001,-73.8951]';
% siteTable = table(names,devices,lats,lons);
cmap = lines(length(unique(devices)));

figure(99),clf
h = geoscatter(lats,lons,200,cmap(devices,:),'filled')
geobasemap topographic
geolimits([41.01,41.06],[-73.95,-73.84]);
% hold on % workaround to display discrete color legend instead of colorbar
% for i=1:size(cmap,1)
% h(i)=geoscatter(lats(1),lons(1),0.5,cmap(i,:),'filled');
% end
% legend(h, unique(devices))
% hold off
t=h.Parent;
t.LatitudeLabel.String="";
t.LatitudeAxis.FontSize=25;
t.LongitudeLabel.String="";
t.LongitudeAxis.FontSize=25;
t.LongitudeAxis.TickLabelFormat = '-dd';
t.LatitudeAxis.TickLabelFormat = '-dd';

savename = 'P:\users\cohen_rebecca_rec297\CCB\HudsonProject\PiermontSites.png';
saveas(gcf,savename,'png');


%% Plot shortnose sturgeon monitoring site and colocated telemetry receiver

names = {'Shortnose','550786','GillNet','GillNet','GillNet'};
device = categorical({'ST300','Receiver 550786','Gill Net','Gill Net','Gill Net'});
% lats = [41.86768,41.86891036,41.8679653684949,41.8671570351114,41.86636537];
% lons = [-73.93382,-73.93388865,-73.9350936537254,-73.93504865,-73.93500032];
lats = [41.86768,41.86891036,41.8762903839038,41.87726205,41.87797539];
lons = [-73.93382,-73.93388865,-73.93894032,-73.93923699,-73.93950365 ];
siteTable = table(names,device,lats,lons);

figure(999)
gb1 = geobubble(siteTable.lats,siteTable.lons,1,siteTable.device);
gb1.BubbleWidthRange = 10;
gb1.Basemap = 'topographic';
gb1.ColorLegendTitle = 'Device';
geolimits([41.85,41.88],[-73.9750 -73.9083])



%% Plot data coverage
siteNames = {'Norrie Point','NP';'Atlantic Sturgeon Spawning','AS';'Black Creek','BC';'Enderkill Creek','EC';...
    'Furnacebrook Creek, Above Dam','FBAD';'Furnacebrook Creek, Below Dam','FBBD';'Shortnose Sturgeon Overwintering','SS';...
    'Tivoli Bay North','TBN';'Tivoli Bay South','TBS';'Tivoli Channel','TCh';'Iona Island','II';'Iona Channel','ICh';...
    'Piermont Marsh','PM';'Piermont Channel','PCh';'Stockport Creek','SCr';'Stockport Channel','SCh'};

deployInfo = readtable('W:\projects\2022_NOAA-NERRS_HudsonNY_144488\DeploymentInfo.xlsx','VariableNamingRule','preserve');
for i=1:height(deployInfo)
    if strncmp(deployInfo.Site(i),'',10)
        deployInfo.Site(i) = deployInfo.Site(i-1);
    end
end
deployInfo([28,54,83:84],:) = [];
% sites = unique(deployInfo.Site);
siteLats = [];
siteLons = [];
for i=1:size(siteNames,1)
    siteInds = find(contains(deployInfo.Deployment,siteNames{i,2}));
    siteLats(i) = mean(deployInfo.Lat(siteInds),'omitnan');
    siteLons(i) = mean(deployInfo.Lon(siteInds),'omitnan');
end

% Sort sites south->north
[sortedLats, I] = sort(siteLats);
sortedSites = siteNames(I,:);

% Convert deployment start/end dates to plottable format
startDates = datetime(deployInfo.("Start Datetime"),'InputFormat','dd-MM-yyyy HH:mm:ss');
stopDates = datetime(deployInfo.("End Datetime"),'InputFormat','dd-MM-yyyy HH:mm:ss');

figure(77),clf
cmap = hsv(size(sortedSites,1));
hold on
for i=1:length(sortedSites)
    siteInds = find(contains(deployInfo.Deployment,sortedSites{i,2}));
    dataRanges = [(startDates(siteInds))';(stopDates(siteInds))'];
    yval = (ones(length(siteInds),1)*i)';
%     plot(dataRanges,[yval;yval],'-','Color',cmap(i,:),'LineWidth',10);
    plot(dataRanges,[yval;yval],'-','Color',[0,102,204]./255,'LineWidth',10);
end
hold off
datetick('x',12)
ylim([0.5,size(sortedSites,1)+0.5])
xlim([min(startDates),max(stopDates)])
yticks([1:size(sortedSites,1)]);
yticklabels(sortedSites(:,1));
% ytickangle(45)
set(gca,'FontSize',20)

savename = 'P:\users\cohen_rebecca_rec297\CCB\HudsonProject\AllHudson_DataCoverage_Plain.png';
saveas(gcf,savename,'png');


%% Plot just NERRS data coverage
siteNames = {'Atlantic Sturgeon Spawning','AS';'Black Creek','BC';...
    'Shortnose Sturgeon Overwintering','SS';'Tivoli Bay North','TBN';'Tivoli Bay South','TBS';...
    'Tivoli Channel','TCh';'Iona Island','II';'Iona Channel','ICh';'Piermont Marsh','PM';...
    'Piermont Channel','PCh';'Stockport Creek','SCr';'Stockport Channel','SCh'};

deployInfo = readtable('W:\projects\2022_NOAA-NERRS_HudsonNY_144488\DeploymentInfo.xlsx','VariableNamingRule','preserve');
for i=1:height(deployInfo)
    if strncmp(deployInfo.Site(i),'',10)
        deployInfo.Site(i) = deployInfo.Site(i-1);
    end
end
deployInfo([39,65],:) = [];
% sites = unique(deployInfo.Site);
siteLats = [];
siteLons = [];
for i=1:size(siteNames,1)
    siteInds = find(contains(deployInfo.Deployment,siteNames{i,2}));
    siteLats(i) = mean(deployInfo.Lat(siteInds),'omitnan');
    siteLons(i) = mean(deployInfo.Lon(siteInds),'omitnan');
end

% Sort sites south->north
[sortedLats, I] = sort(siteLats);
sortedSites = siteNames(I,:);

% Convert deployment start/end dates to plottable format
startDates = datetime(deployInfo.("Start Datetime"),'InputFormat','dd-MM-yyyy HH:mm:ss');
stopDates = datetime(deployInfo.("End Datetime"),'InputFormat','dd-MM-yyyy HH:mm:ss');

figure(77),clf
cmap = hsv(size(sortedSites,1));
hold on
for i=1:length(sortedSites)
    siteInds = find(contains(deployInfo.Deployment,sortedSites{i,2}));
    dataRanges = [(startDates(siteInds))';(stopDates(siteInds))'];
    yval = (ones(length(siteInds),1)*i)';
%     plot(dataRanges,[yval;yval],'-','Color',cmap(i,:),'LineWidth',10);
    plot(dataRanges,[yval;yval],'-','Color',[0,102,204]./255,'LineWidth',10);
end
hold off
ylim([0.5,size(sortedSites,1)+0.5])
% xticks([datetime(2023,1,1),datetime(2023,4,1),datetime(2023,7,1),datetime(2023,10,1)])
xlim([datetime(2023,1,1,0,0,0),datetime(2024,10,1,0,0,0)])
yticks([1:size(sortedSites,1)]);
yticklabels(sortedSites(:,1));
% ytickangle(45)
set(gca,'FontSize',20)

savename = 'P:\users\cohen_rebecca_rec297\CCB\HudsonProject\NERRS_DataCoverage.png';
saveas(gcf,savename,'png');


%% Plot Atlantic Sturgeon coverage
siteNames = {'Atlantic Sturgeon Spawning','AS'};
deployInfo = readtable('W:\projects\2022_NOAA-NERRS_HudsonNY_144488\DeploymentInfo.xlsx','VariableNamingRule','preserve');
for i=1:height(deployInfo)
    if strncmp(deployInfo.Site(i),'',10)
        deployInfo.Site(i) = deployInfo.Site(i-1);
    end
end
deployInfo([28,54,83:84],:) = [];
% sites = unique(deployInfo.Site);
siteLats = [];
siteLons = [];
for i=1:size(siteNames,1)
    siteInds = find(contains(deployInfo.Deployment,siteNames{i,2}));
    siteLats(i) = mean(deployInfo.Lat(siteInds),'omitnan');
    siteLons(i) = mean(deployInfo.Lon(siteInds),'omitnan');
end

% Sort sites south->north
[sortedLats, I] = sort(siteLats);
sortedSites = siteNames(I,:);

% Convert deployment start/end dates to plottable format
startDates = datetime(deployInfo.("Start Datetime"),'InputFormat','dd-MM-yyyy HH:mm:ss');
stopDates = datetime(deployInfo.("End Datetime"),'InputFormat','dd-MM-yyyy HH:mm:ss');

figure(77),clf
cmap = hsv(size(sortedSites,1));
hold on
for i=1:size(sortedSites,1)
    siteInds = find(contains(deployInfo.Deployment,sortedSites{i,2}));
    dataRanges = [(startDates(siteInds))';(stopDates(siteInds))'];
    yval = (ones(length(siteInds),1)*i)';
%     plot(dataRanges,[yval;yval],'-','Color',cmap(i,:),'LineWidth',10);
    plot(dataRanges,[yval;yval],'-','Color',[0,102,204]./255,'LineWidth',10);
end
hold off
xlim([datetime(2021,6,1,0,0,0),datetime(2021,7,20,0,0,0)])
yticks([])
title('Data Coverage')
set(gca,'FontSize',20)

savename = 'P:\users\cohen_rebecca_rec297\CCB\HudsonProject\NERRS_DataCoverage.png';
saveas(gcf,savename,'png');

