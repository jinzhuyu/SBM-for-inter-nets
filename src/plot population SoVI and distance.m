% run after "plot the probability of interlinks to water nodes and power
% nodes
waterNodePopuNorm = WaterPowerNodeTable.PopulNorm(bitwise_equal(WaterPowerNodeTable.NodeLabel,waterLabelSort)); % return bitwise index

figure();
bar(waterNodePopuNorm,'FaceColor',[0.1 0.05 0.6]) 
% dark blue [0.1 0.05 0.6]; dark red [126 0 6]/256
%c = uisetcolor %http://www.thebrooklynpress.com/colors
ylim([0 1])
xticks(1:length(waterNodePopuNorm))
% bar(nodes_can_Fail,nodeImportance)
set(gca,'XTickLabel',waterLabelSort)%,'YTick',0:15)

% xt = get(gca, 'XTick');
set(gca, 'FontSize', fontsize_tick)

xlabel ('Water node number','fontweight','bold','fontsize',fontsize_label)
ylabel ('Population (normalized)','fontweight','bold','fontsize',fontsize_label)

x0=500+i*10;
y0=300+i*10;
wthFig=1200;
htFig=wthFig*0.25;
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



% population at power nodes
PowerNodePopuNorm = WaterPowerNodeTable.PopulNorm(bitwise_equal(WaterPowerNodeTable.NodeLabel,powerLabelSort));

figure();
bar(PowerNodePopuNorm,'FaceColor',[126 0 6]/256) 
% dark blue [0.1 0.05 0.6]; dark red [126 0 6]/256
%c = uisetcolor %http://www.thebrooklynpress.com/colors
ylim([0 1])
xticks(1:length(PowerNodePopuNorm))
% bar(nodes_can_Fail,nodeImportance)
set(gca,'XTickLabel',powerLabelSort-100)%,'YTick',0:15)

% xt = get(gca, 'XTick');
set(gca, 'FontSize', fontsize_tick)

xlabel ('Power node number','fontweight','bold','fontsize',fontsize_label)
ylabel ('Population (normalized)','fontweight','bold','fontsize',fontsize_label)

x0=500+i*10;
y0=300+i*10;
wthFig=1200;
htFig=wthFig*0.25;
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


% SoVI
% water SoVI
waterNodeSovi = WaterPowerNodeTable.Sovi(bitwise_equal(WaterPowerNodeTable.NodeLabel,waterLabelSort));

figure();
bar(waterNodeSovi,'FaceColor',[0.1 0.05 0.6]) 
% dark blue [0.1 0.05 0.6]; dark red [126 0 6]/256
%c = uisetcolor %http://www.thebrooklynpress.com/colors
ylim([0 1])
xticks(1:length(waterNodeSovi))
% bar(nodes_can_Fail,nodeImportance)
set(gca,'XTickLabel',waterLabelSort)%,'YTick',0:15)

% xt = get(gca, 'XTick');
set(gca, 'FontSize', fontsize_tick)

xlabel ('Water node number','fontweight','bold','fontsize',fontsize_label)
ylabel ('SoVI','fontweight','bold','fontsize',fontsize_label)

x0=500+i*10;
y0=300+i*10;
wthFig=1200;
htFig=wthFig*0.25;
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



% SoVI at power nodes
PowerNodeSovi = WaterPowerNodeTable.Sovi(bitwise_equal(WaterPowerNodeTable.NodeLabel,powerLabelSort));

figure();
bar(PowerNodeSovi,'FaceColor',[126 0 6]/256) 
% dark blue [0.1 0.05 0.6]; dark red [126 0 6]/256
%c = uisetcolor %http://www.thebrooklynpress.com/colors
ylim([0 1])
xticks(1:length(PowerNodeSovi))
% bar(nodes_can_Fail,nodeImportance)
set(gca,'XTickLabel',powerLabelSort-100)%,'YTick',0:15)

% xt = get(gca, 'XTick');
set(gca, 'FontSize', fontsize_tick)

xlabel ('Power node number','fontweight','bold','fontsize',fontsize_label)
ylabel ('SoVI','fontweight','bold','fontsize',fontsize_label)

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% distance
% normalized the reciprical of distance water
distRecipFromWaterNodesKm = 1./distFromWaterNodesKm;
distRecipFromWaterNodesKmNorm = truncated_normalize(distRecipFromWaterNodesKm,normLowBound);
distRecipFromWaterNodesKmNormRowNorm = bsxfun(@rdivide,distRecipFromWaterNodesKmNorm,...
                                              sum(distRecipFromWaterNodesKmNorm,2));
distRecipFromWaterNodesKmNormMean = mean(distRecipFromWaterNodesKmNormRowNorm);
% % plot distribution of pr of water nodes to 9th gate station
% bar(distRecipFromWaterNodesKmNorm(9,:)./sum(distRecipFromWaterNodesKmNorm(9,:)))
% max(distRecipFromWaterNodesKmNorm(9,:)./sum(distRecipFromWaterNodesKmNorm(9,:))) % 0.3610
% water nodes

waterNodeDistRecipNorm = distRecipFromWaterNodesKmNormMean(...
                            bitwise_equal(WaterNodeInterDelivery,waterLabelSort)); % return bitwise index

figure();
temp_water_dist = sort(truncated_normalize(waterNodeDistRecipNorm,0.01),'descend');
bar(temp_water_dist,'FaceColor',[0.1 0.05 0.6]) 
% dark blue [0.1 0.05 0.6]; dark red [126 0 6]/256
%c = uisetcolor %http://www.thebrooklynpress.com/colors
ylim([0 1])
xticks(1:length(waterNodeDistRecipNorm))
% bar(nodes_can_Fail,nodeImportance)
waterLabelSort_temp = [23,21,49,29,34,30,27,25,39,42,31,33,35,37,18,28,38,...
                  22,19,41,16,40,26,36,44,20,45,48,46,47,32,43,24,17];
set(gca,'XTickLabel',waterLabelSort_temp)%,'YTick',0:15)

% xt = get(gca, 'XTick');
set(gca, 'FontSize', fontsize_tick)

xlabel ('Water node number','fontweight','bold','fontsize',fontsize_label)
ylabel ({'Reciprocal of distance';'(normalized)'},'fontweight','bold','fontsize',fontsize_label)

x0=500+i*10;
y0=300+i*10;
wthFig=1200;
htFig=wthFig*0.25;
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


% power
% normalized the reciprical of distance

distRecipFromPowerNodesKm = 1./distFromPowerNodesKm;
distRecipFromPowerNodesKmNorm = truncated_normalize(distRecipFromPowerNodesKm,normLowBound);
distRecipFromPowerNodesKmNormRowNorm = bsxfun(@rdivide,distRecipFromPowerNodesKmNorm,...
                                                sum(distRecipFromPowerNodesKmNorm,2));                         
distRecipFromPowerNodesKmNormMean = mean(distRecipFromPowerNodesKmNormRowNorm); % from each PowerNodeInterDelivery 16-49

powerNodedistRecipNorm = distRecipFromPowerNodesKmNormMean(...
                         bitwise_equal(PowerNode12Substation,powerLabelSort)); % return bitwise index

figure();

temp_power_dist = sort(truncated_normalize(powerNodedistRecipNorm,0.01),'descend');
bar(temp_power_dist,'FaceColor',[126 0 6]/256)
% dark blue [0.1 0.05 0.6]; dark red [126 0 6]/256
%c = uisetcolor %http://www.thebrooklynpress.com/colors
ylim([0 1])
xticks(1:length(powerNodedistRecipNorm))
% bar(nodes_can_Fail,nodeImportance)
powerLabelSort_temp = [32,37,27,33,28,35,31,30,34,26,41,36,39,20,43,40,44,42,45,38]+100
set(gca,'XTickLabel',powerLabelSort_temp-100)%,'YTick',0:15)

% xt = get(gca, 'XTick');
set(gca, 'FontSize', fontsize_tick)

xlabel ('Power node number','fontweight','bold','fontsize',fontsize_label)
ylabel ({'Reciprocal of distance';'(normalized)'},'fontweight','bold','fontsize',fontsize_label)

x0=500+i*10;
y0=300+i*10;
wthFig=1200;
htFig=wthFig*0.25;
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


