% distribution of resilience at some time point
lWidth=3;
fontsize_label=14;
fontsize_tick=12;

t_sel=[5 30 50];

for ii=1:length(t_sel)
%     title=sprintf('PDF of the Resilience at %d',t_sel(ii));
%     PDF_resil_stat=figure('name',title);
    figure()
    hold on;
    hist(ResilAll(:,t_sel(ii)))
    h1 = findobj(gca,'Type','patch');
    set(h1,'facealpha',0.85,'FaceColor',[1 1 1].*150/255);
%     h1.FaceColor=[1 1 1].*150/255;

    hist(ResilTotDynamicAll(:,t_sel(ii)));
    h2 = findobj(gca,'Type','patch');
%     h2.FaceColor=[0 0 1];
    set(h2,'facealpha',0.85)
    hold off;
    
%   [pdf_resil_stat,pdf_resil_stat_xi] = ksdensity(ResilAll(:,t_sel(ii)));
%   trapz(pdf_meanService_xi,pdf_meanService)=0.9999
%   plot(pdf_resil_stat_xi,pdf_resil_stat/sum(pdf_resil_stat),'LineWidth',lWidth,'Color','b')

    
%     hold on;
%     [pdf_resil_dyn,pdf_resil_dyn_xi] = ksdensity(ResilTotDynamicAll(:,t_sel(ii)));
%     %     trapz(pdf_meanService_xi,pdf_meanService)=0.9999
%     plot(pdf_resil_dyn_xi,pdf_resil_dyn/sum(pdf_resil_dyn),'LineWidth',lWidth,'Color','k')
%     hold off;
    
    xlabel_text=sprintf('Resilience at time step %d',t_sel(ii));
    xlabel(xlabel_text,'FontSize',fontsize_label,'FontWeight','bold')
    ylabel('Frequency','FontSize',fontsize_label,'FontWeight','bold')
    legend('Static ranking','Dynamic ranking','location','northwest');
    
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', fontsize_tick)
    % set(gcf, 'Color', 'None') % make the background transparent
%     grid on;
%     grid minor
    ax = gca;
    ax.YGrid = 'on';
    ax.GridLineStyle = '-';

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


