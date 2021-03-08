Water2PowerInterlinkSample_water_nodes_3d = Water2PowerInterlinkSample(:,1,:);
Water2PowerInterlinkSample_water_nodes_2d = squeeze(Water2PowerInterlinkSample_water_nodes_3d);
Water2PowerInterlinkSample_water_nodes = Water2PowerInterlinkSample_water_nodes_2d(:);

% plot node importance
% water node
fontsize_label = 14;
fontsize_tick = 12.5;

% importance by change in the mean
figure();
[freq_water,water_nodes] = hist(Water2PowerInterlinkSample_water_nodes,unique(Water2PowerInterlinkSample_water_nodes));

[freq_water_sort,water_nodes_Index_sort] = sort(freq_water,'descend');
waterLabelSort = water_nodes(water_nodes_Index_sort);
bar(freq_water_sort/(nSampleInterlink*nWaterNodePumpStation)*100,'FaceColor',[0.1 0.05 0.6]) 
% dark blue [0.1 0.05 0.6]; dark red [126 0 6]/256
%c = uisetcolor %http://www.thebrooklynpress.com/colors
ylim([0 8])
xticks(1:length(waterLabelSort))
% bar(nodes_can_Fail,nodeImportance)
set(gca,'XTickLabel',waterLabelSort)%,'YTick',0:15)

% xt  =  get(gca, 'XTick');
set(gca, 'FontSize', fontsize_tick)

xlabel ('Water node number','fontweight','bold','fontsize',fontsize_label)
ylabel ('Probability (×10^{-2})','fontweight','bold','fontsize',fontsize_label)

x0 = 500+i*10;
y0 = 300+i*10;
wthFig = 1200;
htFig = wthFig*0.25;
set(gcf,'position',[x0,y0,wthFig,htFig])

% grid on;
% grid minor;
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';
box off;

ax2 = axes('Position',get(gca,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k');
set(ax2,'YTick', []);
set(ax2,'XTick', []);
box on;



% power nodes
Power2WaterInterlinkSample_water_nodes_3d = Power2WaterInterlinkSample(:,1,:);
Power2WaterInterlinkSample_water_nodes_2d = squeeze(Power2WaterInterlinkSample_water_nodes_3d);
Power2WaterInterlinkSample_water_nodes = Power2WaterInterlinkSample_water_nodes_2d(:);

figure();
[freq_power,power_nodes] = hist(Power2WaterInterlinkSample_water_nodes,unique(Power2WaterInterlinkSample_water_nodes));
% importance by change in the mean
[freq_power_sort,power_nodes_Index_sort] = sort(freq_power,'descend');
powerLabelSort = power_nodes(power_nodes_Index_sort);
bar(freq_power_sort/(nSampleInterlink*nWaterNodePumpStation)*100,'FaceColor',[126 0 6]/256) 
% dark blue [0.1 0.05 0.6]; dark red [126 0 6]/256
% c = uisetcolor %http://www.thebrooklynpress.com/colors
ylim([0 15])
xticks(1:length(powerLabelSort))
% bar(nodes_can_Fail,nodeImportance)
set(gca,'XTickLabel', powerLabelSort-100)%,'YTick',0:15)

% xt = get(gca, 'XTick');
set(gca, 'FontSize', fontsize_tick)

xlabel ('Power node number','fontweight','bold','fontsize',fontsize_label)
ylabel ('Probability (×10^{-2})','fontweight','bold','fontsize',fontsize_label)

x0 = 500+i*10;
y0 = 300+i*10;
wthFig = 1200;
htFig = wthFig*0.25;
set(gcf,'position',[x0,y0,wthFig,htFig])

% grid on;
% grid minor;
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';
box off;

ax2 = axes('Position',get(gca,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k');
set(ax2,'YTick', []);
set(ax2,'XTick', []);
box on;