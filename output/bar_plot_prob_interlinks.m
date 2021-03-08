% bar plot of probability of interlinks
for i=1:4
    prob=PrPower2Water(i,:); %PrWater2Power(i,:);
    nodeLabels=PowerNode12Substation-100; %WaterNodeInterDelivery;

    % plot node importance
    fontsize_label=10;
    fontsize_tick=10;

    % importance by change in the mean
    figure();
    [ProbSort,ProbIndexSort]=sort(prob,'descend');
    labels_sort=nodeLabels(ProbIndexSort);
    bar(ProbSort.*100,'FaceColor',[126 0 6]/256) 
    % dark blue [0.1 0.05 0.6]; dark red [126 0 6]/256
    %c = uisetcolor %http://www.thebrooklynpress.com/colors
    ylim([0 0.18]*100)
    xticks(1:length(labels_sort))
    % bar(nodes_can_Fail,nodeImportance)
    set(gca,'XTickLabel',labels_sort)%,'YTick',0:15)

    % xt = get(gca, 'XTick');
    set(gca, 'FontSize', fontsize_tick)

    xlabel ('Power node number','fontweight','bold','fontsize',fontsize_label)
    ylabel ('Probability (×10^{-2})','fontweight','bold','fontsize',fontsize_label)

    x0=500+i*10;
    y0=300+i*10;
    wthFig=700;
    htFig=wthFig*1.85/8;
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
end
