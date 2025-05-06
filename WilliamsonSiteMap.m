% Plot recording sites colored by project

projects = categorical({'BOEM','NYSDEC','NEAq'}); % name of project each device belongs to
lats = [41.1417465,40.3477982,38.30304]'; % lat of each device
lons = [-71.1041718,-71.224167,-74.65363]'; % lon of each device
saveDir = ''; % path to folder where you'd like to save the map
saveName = 'ThesisSites'; % what to name the saved map

%%
cmap = lines(length(unique(devices)));

figure(99),clf
geobasemap landcover % set basemap
h = geoscatter(lats,lons,200,cmap(projects,:),'filled'); % plot recording locations colored by project
geolimits([34.5,43.5],[-77,-68]) % set lat/lon bounds of map (also depends on size of figure window)
hold on % workaround to display discrete color legend instead of colorbar
for i=1:size(cmap,1)
h(i)=geoscatter(lats(1),lons(1),0.5,cmap(i,:),'filled');
end
legend(h, unique(projects))
hold off
t=h.Parent;
t.LatitudeLabel.String=''; % set y-axis label; set to '' to make the label empty
t.LatitudeAxis.FontSize=10; % set font size of y-axis tick labels
t.LongitudeLabel.String=''; % set x-axis label; set to '' to make the label empty
t.LongitudeAxis.FontSize=10; % set font size of x-axis tick labels
t.LongitudeAxis.TickLabelFormat = '-dd'; % set format of lat/lon values
t.LatitudeAxis.TickLabelFormat = '-dd';% set format of lat/lon values

saveas(gcf,fullfile(saveDir,saveName),'png');
exportgraphics(gcf,fullfile(saveDir,saveName),'ContentType','vector');

