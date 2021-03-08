function [h1,h2]=plot_with_bounds(x,y,colorFill,colorMean,colorBounds)
% Plot mean value with lower and upper bounds
% input: 
% x: 1d array; y: 2d array;

% https://www.mathworks.com/matlabcentral/answers/180829-shade-area-between-graphs
% no fill between bounds
yMean = mean(y);
% yVar = var(y);
% z=1.96;
% sigma=std(y);
% yLow=yMean-z*sigma/sqrt(length(y(:,1)));
% yUp=yMean+z*sigma/sqrt(length(y(:,1)));
% yLow = yMean - 2*yVar;
% yUp = yMean + 2*yVar;
% yMin= min(y);
% find the second smallest element of each column in y
y_copy = y;

% % repace zeros in each column with the mean of each column
% for i=1:2
%     for iCol=1:length(y_copy(1,:))
%         y_copy_icol = y_copy(:,iCol);
%     %     if min(y_copy_icol)==0
%         y_copy_icol(y_copy_icol==min(y_copy_icol)) = mean(y_copy(:,iCol));
%     %     else
%     %         y_copy_icol(y_copy_icol==min(y_copy_icol)) = mean(y_copy(:,iCol));        
%     %     end
%         y_copy(:,iCol) = y_copy_icol;
%     end
% end
% for i=1:6
%     y_copy(bsxfun(@eq, y_copy, min(y_copy))) = y_mean; % element wise binary operation. replace min with max
% end

% lower and upper bounds of y
% yLow = min(y_copy);
% yUp = max(y);
% if

% yLow = quantile(y_copy,0.00001);
yUp = max(y_copy);

yLow = min(y_copy);
for ii = 1:43

  yLow(ii) = quantile(y_copy(:,ii),0.0125);
end
% % x0=500;
% y0=300;
% wthFig=550;
% htFig=wthFig*5.5/8;
% set(gcf,'position',[x0,y0,wthFig,htFig])
% % figure()
% hold on
% scatter(x, y, 'k')
% fill between bounds
xFill=[x,fliplr(x)];
yFill=[yLow,fliplr(yUp)];

h2 = fill(xFill,yFill,colorFill);  %'EdgeColor','none'
set(h2,'FaceAlpha',0.2)
% Plot yMean
h1 = plot(x, yMean, 'Color', colorMean, 'LineStyle', '-', 'LineWidth', 1.75);

widthBnd =1.25;
% sizeLabel=14;
% sizeTick=12;
% % Plot upper and lower bounds, calculated as 0.3 from yfit
h3=plot(x, yUp, 'Color', colorBounds, 'LineStyle', '-', 'LineWidth', widthBnd-0.15);
% % alpha(h2,.1)
h4=plot(x, yLow, 'Color', colorBounds, 'LineStyle', '-', 'LineWidth', widthBnd);
% % alpha(h3,.1)
% % set(get(get(h4(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% % Add a legend and axis labels
% hLegend=legend([h1 h4],'Mean', 'Range', 'Location', 'Southeast');
% set(hLegend,'FontSize',sizeTick)
% 
% xlabel('Time step','FontSize',sizeLabel,'fontweight','bold')
% ylabel('Resilience','FontSize',sizeLabel,'fontweight','bold')
% xlim([0 70])
% ylim([0 1.05])
% hold off
% 
% % xt = get(gca, 'XTick');
% % set(gca, 'FontSize', sizeTick)
% % % set(gcf, 'Color', 'None') % make the background transparent
% % grid on;
% % %     grid minor
% % box off;
% % ax2 = axes('Position',get(gca,'Position'),...
% %            'XAxisLocation','top',...
% %            'YAxisLocation','right',...
% %            'Color','none',...
% %            'XColor','k','YColor','k');
% % set(ax2,'YTick', []);
% % set(ax2,'XTick', []);
% % box on;
% % hold off
end