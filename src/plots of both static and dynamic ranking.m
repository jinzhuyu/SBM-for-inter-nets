% plot the resilience curve
    fontsize_label=12;
    fontsize_tick=11;
    lWidth=1.0;
    lWidth2=1.75;
    
    figure
    hold on
    plot(1:nSteps,ResilWater,':k+','LineWidth',lWidth)
    plot(1:nSteps,Resil,'-c*','LineWidth',lWidth)
    plot(1:nSteps,ResilPower,'--mx','LineWidth',lWidth)
    
    plot(1:nSteps,ResilWaterDynamic,'-.b','LineWidth',lWidth2)
    plot(1:nSteps,ResilTotDynamic,'-g','LineWidth',lWidth2)
    plot(1:nSteps,ResilPowerDynamic,'--r','LineWidth',lWidth2)
%     ylim([min(Resil)-0.05 1.05])
    
    xlabel('Time step','FontSize',fontsize_label,'FontWeight','bold')
    ylabel('Resilience','FontSize',fontsize_label,'FontWeight','bold')
    legend({'Water network based on static ranking','Mean value based on static ranking','Power network based on static ranking',...
        'Water network based on dynamic ranking','Mean value based on dynamic ranking','Power network based on dynamic ranking'
        },'Location','southeast','FontSize',fontsize_label)
    
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', fontsize_tick)
    % set(gcf, 'Color', 'None') % make the background transparent
    grid on;
%     grid minor
    box off;
    ax2 = axes('Position',get(gca,'Position'),...
               'XAxisLocation','top',...
               'YAxisLocation','right',...
               'Color','none',...
               'XColor','k','YColor','k');
    set(ax2,'YTick', []);
    set(ax2,'XTick', []);
    box on;
    hold off