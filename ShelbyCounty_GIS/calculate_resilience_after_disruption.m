function ResiTotInitialTemp=calculate_resilience_after_disruption(WaterPowerNodeTableTemp, WaterPowerLinkArray,nScen)

% WaterPowerLinkTable=PowerAndWaterAllLinkSample(:,:,iSample);
% node table: NodeLabel, NodeType, Lat, Long.
% link table: link id, start and end node id, coordinates.
epicenter = [35.3 -90.3];
Mw = 7.7;

DistInKm= distdim(distance(epicenter(1), epicenter(2), WaterPowerNodeTableTemp.Lat, WaterPowerNodeTableTemp.Long), 'deg', 'km');
    
term1 = 3.79;
term2 = 0.298*(Mw - 6);
term3 = 0.0536*(Mw - 6)^2;
term4 = log10(DistInKm);
term5 = 0.00135*DistInKm;
    
PGA = 10.^(term1 + term2 - term3 -term4 - term5); % pga

% waterPowerNodeTable=array2table(waterPowerNode,'VariableNames',{'NodeID','NodeType','NodeIDLabel','Lat','Long','distToCenter','PGA'});

% number of nodes, only the first 15 nodes can fail, though.
% num_of_nodes = waterPowerNodeTable.NodeIDLabel < 3;
NodeIdCanFail = WaterPowerNodeTableTemp.NodeLabel<104;

GravityAcceleration=9.80665;
PGAInG=PGA/GravityAcceleration/100;  
% lambda_PfByPGA = [log(1.5) log(1.5) log(3.0) log(0.47) log(0.7) log(0.9)];   % Table 1, T. Adachi, B.R. Ellingwood / Reliability Engineering and System Safety 93 (2008), page 84     
lambdaPGA = [log(1.5) log(2.0) log(2.5) log(1.2) log(1.3) log(1.4) log(5)];  
% assumed to generate more reasonable; pumping stations have power backup 
zetaPGA = [0.8 0.6 0.2 0.4 0.4 0.4 0.1]; 

WaterPowerNodeTableTemp.NodeType(WaterPowerNodeTableTemp.NodeType>100)=WaterPowerNodeTableTemp.NodeType(WaterPowerNodeTableTemp.NodeType>100)-97;

for iNodeCanFail = 1:size(NodeIdCanFail,1)
    iNodeType = WaterPowerNodeTableTemp.NodeType(iNodeCanFail);
    NodePf=normcdf(PGAInG(iNodeCanFail),lambdaPGA(iNodeType),zetaPGA(iNodeType)); 
    WaterPowerNodeTableTemp.NodePf(iNodeCanFail) = NodePf;    
end

 WaterPowerNodeTableTemp.NodePf(WaterPowerNodeTableTemp.NodeType==7)=0;

% calcualte PGV
WaterPowerLinkTable=array2table(WaterPowerLinkArray,'VariableNames',...
    {'NodeStartId','NodeEndId','TotLength','StartNodeLat',...
    'StartNodeLong','EndNodeLat','EndNodeLong'});


epicenter = [35.3 -90.3];
Mw = 7.7;

StartNodeDist_km = distdim(distance(epicenter(1), epicenter(2),WaterPowerLinkTable.StartNodeLat,...
    WaterPowerLinkTable.StartNodeLong),'deg','km');
EndNodeDist_km = distdim(distance(epicenter(1), epicenter(2), WaterPowerLinkTable.EndNodeLat, ...
    WaterPowerLinkTable.EndNodeLong),'deg','km');
     
term1 = 2.04;
term2 = (0.422*(Mw - 6));
term3 = (0.0373*(Mw - 6)^2);

StartNodeTerm4 = log10(StartNodeDist_km);
EndNodeTerm4 = log10(EndNodeDist_km);
    
WaterPowerLinkTable.StartNodePGV_ins= 10.^(term1 + term2 - term3 -StartNodeTerm4)./2.54; 
WaterPowerLinkTable.EndNodePGV_ins = 10.^(term1 + term2 - term3 -EndNodeTerm4)./2.54;     

% links
aLink=[0.002 0.001]; % 0.00187;
KLink=0.5;

WaterPowerLinkTable.NodeStartType=zeros(length(WaterPowerLinkTable.NodeStartId),1);
WaterPowerLinkTable.NodeStartType(WaterPowerLinkTable.NodeStartId<100)=1; % water link
WaterPowerLinkTable.NodeStartType(WaterPowerLinkTable.NodeStartId>100)=2;
for iLink = 1:length(WaterPowerLinkTable.NodeStartId)

    RRWaterLink=aLink(WaterPowerLinkTable.NodeStartType(iLink)).*KLink.*...
                     (WaterPowerLinkTable.StartNodePGV_ins(iLink)+...
                     WaterPowerLinkTable.EndNodePGV_ins(iLink))./2;

    WaterPowerLinkTable.LinkPf(iLink) = 1-exp(-RRWaterLink.*WaterPowerLinkTable.TotLength(iLink).*3280.84./1000);

end

%----------------------------
% Evaluate network resilience

nodeTable = WaterPowerNodeTableTemp;
linkTable =  WaterPowerLinkTable;

nodeTable.NodeInterdependType=strings(length(nodeTable.NodeLabel),1);

linkTable.LinkType=strings(length(linkTable.NodeStartId),1);

linkTable.LinkType(linkTable.NodeStartId<100 & linkTable.NodeEndId<100)="WW";
linkTable.LinkType(linkTable.NodeStartId>100 & linkTable.NodeEndId>100)="PP";
linkTable.LinkType(linkTable.NodeStartId<100 & linkTable.NodeEndId>100)="WP";
linkTable.LinkType(linkTable.NodeStartId>100 & linkTable.NodeEndId<100)="PW";

nodeTable.NodeInterdependType(ismember(nodeTable.NodeLabel,linkTable.NodeStartId(linkTable.LinkType=='WP')))="WaterSupplyPower";
nodeTable.NodeInterdependType(ismember(nodeTable.NodeLabel,linkTable.NodeStartId(linkTable.LinkType=='PW')))="PowerSupplyWater";
nodeTable.NodeInterdependType(ismember(nodeTable.NodeLabel,linkTable.NodeEndId(linkTable.LinkType=='WP')))="PowerDemandWater";
nodeTable.NodeInterdependType(ismember(nodeTable.NodeLabel,linkTable.NodeEndId(linkTable.LinkType=='PW')))="WaterDemandPower";

% create base case u-v pairs: the ID of starting and end nodes respectively
s_base = linkTable.NodeStartId(:);
t_base = linkTable.NodeEndId(:);

% generate grpah
G = graph();
G = addnode(G, size(nodeTable.NodeLabel,1));
G = addedge(G, s_base,t_base);
% figure(1)
% plot(G)
% d = distances(G);

% generate u vectors, size of edges + links
% double loop Monte Carlo 

nodes_can_Fail_IDs = nodeTable.NodeType <8;
nodes_can_Fail = nodeTable.NodeLabel(nodes_can_Fail_IDs);

num_nodes = size(nodes_can_Fail,1);
num_links = size(linkTable,1);

U_nodes = rand(nScen, num_nodes);
U_links = rand(nScen, num_links);% generate num_scenarios*num_links vectors of random #

probFailure_nodes = repmat(nodeTable.NodePf(nodes_can_Fail_IDs)',nScen,1);
probFailure_links = repmat(linkTable.LinkPf',nScen,1);

% if U < pf_fail, then the node/link is disabled 
NodesFailedLogic = U_nodes <probFailure_nodes;  % acceptance-rejection sampling to generate pf*num_scenarios
LinksFailedLogic = U_links <probFailure_links;

% Functionality_xi = zeros(numScen, max(nodeTable.NodeLabel));
% nSteps=5;
% TimeSteps= linspace(1,nSteps,nSteps);

ResiTot = zeros(nScen, nSteps);
ResiWater=zeros(nScen, nSteps);
ResiPower=zeros(nScen, nSteps);

global WaterDemandNodeLabel;
WaterDemandNodeLabel = nodeTable.NodeLabel(nodeTable.NodeType == 3);  % 16-49
global PowerDemandNodeLabel;
PowerDemandNodeLabel = nodeTable.NodeLabel(nodeTable.NodeType == 6);  % | nodeTable.NodeType == 5);
global nWaterDemandNode;
global nPowerDemandNode;
nWaterDemandNode=length(WaterDemandNodeLabel);
nPowerDemandNode=length(PowerDemandNodeLabel);
% numNodeIdAll=max(nodeTable.NodeLabel);

LinksOrig = [s_base t_base]; % Original links
NodesCanFailOrigLabel = nodes_can_Fail'; % Original nodes


for iScen = 1:nScen    
    
    % initial resilience
    
    linkTable.FailureState=LinksFailedLogic(iScen,:)';
    nodeTable.FailureState=NodesFailedLogic(iScen,:)';
    
    % find failed power node supply to water
    PowerFailToSupplyWaterNode_Label=nodeTable.NodeLabel(nodeTable.NodeInterdependType...
                                    =="PowerSupplyWater" & nodeTable.FailureState==1);
 
    WaterFailToSupplyPowerNode_Label=nodeTable.NodeLabel(nodeTable.NodeInterdependType...
                                    =="WaterSupplyPower" & nodeTable.FailureState==1);
    
    % find the correpsonding water node that demands power
    
    WaterNodeNeedNewPowerNode_Label=linkTable.NodeEndId(ismember(linkTable.NodeStartId,...
                                    PowerFailToSupplyWaterNode_Label) & linkTable.NodeEndId<100);
    PowerNodeNeedNewWaterNode_Label=linkTable.NodeEndId(ismember(linkTable.NodeStartId,...
                                    WaterFailToSupplyPowerNode_Label) & linkTable.NodeEndId>100);
    
    PowerNodeOperational=nodeTable.NodeLabel(ismember(nodeTable.NodeLabel,...
                                    WaterPowerNodeOperational) & nodeTable.NodeType==6); % 6 end user power node type id=6
    PowerNodeOperationalRowLogic=ismember(nodeTable.NodeLabel,PowerNodeOperational);
    WaterNodeOperational=nodeTable.NodeLabel(ismember(nodeTable.NodeLabel,...
                                    WaterPowerNodeOperational) & nodeTable.NodeType==3); % 6 end user power node type id=3
    WaterNodeOperationalRowLogic=ismember(nodeTable.NodeLabel,WaterNodeOperational);
    
    NewPowerSupplyToWaterNodeLabel=zeros(length(WaterNodeNeedNewPowerNode_Label),1);
    NewWaterSupplyToPowerNodeLabel=zeros(length(PowerNodeNeedNewWaterNode_Label),1);
    % find the neareast alternative, operational (nodes in the final link table???) power supply node
    for iWaterNode=1:length(WaterNodeNeedNewPowerNode_Label)
          
        Dist= distdim(distance(nodeTable.Lat(iWaterNode), nodeTable.Long(iWaterNode), ...
            nodeTable.Lat(PowerNodeOperationalRowLogic),nodeTable.Long(PowerNodeOperationalRowLogic)), 'deg', 'km');
        NewPowerSupplyToWaterNodeLabel(iWaterNode,1)=PowerNodeOperational(Dist==min(Dist));  
     
    end
    for iPowerNode=1:length(PowerNodeNeedNewWaterNode_Label)
          
        Dist= distdim(distance(nodeTable.Lat(iPowerNode), nodeTable.Long(iPowerNode), ...
            nodeTable.Lat(PowerNodeOperationalRowLogic),nodeTable.Long(WaterNodeOperationalRowLogic)), 'deg', 'km');
        NewWaterSupplyToPowerNodeLabel(iPowerNode,1)=WaterNodeOperational(Dist==min(Dist));  
     
    end
    
    InterLinksNew = vertcat([NewWaterSupplyToPowerNodeLabel,PowerNodeNeedNewWaterNode_Label],...
                   [NewPowerSupplyToWaterNodeLabel,WaterNodeNeedNewPowerNode_Label]);
    InterLinksNewTable=array2table(InterLinksNew, 'VariableNames',{'NodeStartId','NodeEndId'});
    InterLinksNewTable.FailureState=ones(length(InterLinksNewTable.NodeStartId),1);
    % New link table after building new interdependency links
    LinkTableWithNewInterLinks.NodeStartId = vertcat(linkTable.NodeStartId,InterLinksNewTable.NodeStartId);
    LinkTableWithNewInterLinks.NodeEndId = vertcat(linkTable.NodeEndId,InterLinksNewTable.NodeEndId);                 
    LinkTableWithNewInterLinks.FailureState = vertcat(linkTable.FailureState,InterLinksNewTable.FailureState); 
    
    
    LinkTableWithNewInterLinks=struct2table(LinkTableWithNewInterLinks);
    
%     nSteps=sum(linkTable.FailureState)+sum(nodeTable.FailureState);
%     TimeSteps = linspace(1,nSteps,nSteps);
    
    
    % component importance, will be used as the activity priority 
    
    % operational links, including new interlinks
    LinksNotFailedFinalWithNewInterLinks.NodeStartId = vertcat(LinksNotFailedFinal(:,1),...
                                                       InterLinksNewTable.NodeStartId);
    LinksNotFailedFinalWithNewInterLinks.NodeEndId = vertcat(LinksNotFailedFinal(:,2),...
                                                        InterLinksNewTable.NodeEndId);                 
    LinksNotFailedFinalWithNewInterLinks=struct2table(LinksNotFailedFinalWithNewInterLinks);
    
    % calculate the resilience after disruption with new inter links
    [ResiTotTemp,ResiWaterTemp,ResiPowerTemp] = calculate_resi_given_links(LinksNotFailedFinalWithNewInterLinks);
    ResiWaterWithNewInterLinks=ResiWaterTemp;
    ResiPowerWithNewInterLinks=ResiPowerTemp;
    ResiTotWithNewInterLinks=ResiTotTemp;
    
    % rank the importance of damaged links first
    
    % identify damaged links: failure state==1
    DamagedLinks.NodeStartId = vertcat(linkTable.NodeStartId(linkTable.FailureState==1),InterLinksNewTable.NodeStartId);
    DamagedLinks.NodeEndId = vertcat(linkTable.NodeEndId(linkTable.FailureState==1),InterLinksNewTable.NodeEndId);        
%     DamagedLinks=struct2array(DamagedLinks);
    DamagedLinksTable=struct2table(DamagedLinks);
    DamagedLinksTable.ComponentType(:)="Link";
    % for each damaged link in damaged link table, switch the failure state. failure state=0
    % and append that link to the link not fail final
    
    % sort the importance of nodes
    % identify damaged nodes
    DamagedNodesLabel.NodeStartId=nodeTable.NodeLabel(nodeTable.FailureState);
    DamagedNodesTable=struct2table(DamagedNodesLabel);
    DamagedNodesTable.ComponentType(:)="Node";
    DamagedNodesArray=table2array(DamagedNodesTable(:,1));
    
    DamagedComponents.NodeStartId=vertcat(DamagedLinksTable.NodeStartId,DamagedNodesTable.NodeStartId);
    DamagedComponents.NodeEndId=vertcat(DamagedLinksTable.NodeEndId,DamagedNodesTable.NodeStartId);
    DamagedComponents.ComponentType=vertcat(DamagedLinksTable.ComponentType,DamagedNodesTable.ComponentType);
    DamagedComponentsTable=struct2table(DamagedComponents);
    DamagedComponentsArray=table2array(DamagedComponentsTable(:,1:2));
    % for each failed nodes
    Importancecomponents=zeros(length(DamagedComponentsTable.NodeStartId),1);
    for iDamagedComponent=1:length(DamagedComponentsTable.NodeStartId)
        LinksNotFailedFinalTemp=LinksNotFailedFinal;
        if DamagedComponentsTable.ComponentType(iDamagedComponent)=="Link"
            LinksNotFailedFinalTemp=vertcat(LinksNotFailedFinalTemp(:,:),DamagedComponentsArray(iDamagedComponent,1:2));
            ResiTotTemp=calculate_resi_given_links(LinksNotFailedFinalTemp);
            Importancecomponents(iDamagedComponent)=(ResiTotTemp-ResiTotWithNewInterLinks)/ResiTotWithNewInterLinks;
        else
            NodesFailedScenILabelTemp = NodesCanFailOrigLabel(nodeTable.FailureState & (NodesCanFailOrigLabel~=DamagedNodesArray(iDamagedComponent-height(DamagedLinksTable)))'); % Failed nodes
            LinksStartNodesFailedLogical = ismember(LinksNotFailedScenI(:,1), NodesFailedScenILabelTemp);
            LinksFailedStartNodes = LinksNotFailedScenI(~LinksStartNodesFailedLogical,:);
            LinksFailedEndNodesLogical = ismember(LinksFailedStartNodes(:,2), NodesFailedScenILabelTemp );
            LinksNotFailedFinalTemp = LinksFailedStartNodes(~LinksFailedEndNodesLogical,:);
            
            ResiTotTemp=calculate_resi_given_links(LinksNotFailedFinalTemp);
            Importancecomponents(iDamagedComponent)=(ResiTotTemp-ResiTotWithNewInterLinks)/ResiTotWithNewInterLinks;
    %  length(LinksNotFailedFinalTemp(:,1))
    % Perhaps in calculating the importance, use resilience without new
    % interlinks.
        end
          
    end
    
    [ImportanceComponentsDescend,ImportanceComponentsRowIndex]=sort(Importancecomponents,'descend');
    DamagedComponentsImportanceSort=[DamagedComponentsArray(ImportanceComponentsRowIndex,1:2) ImportanceComponentsDescend]; 

    nSteps=length(DamagedComponentsImportanceSort(:,1));
    LinksNotFailedFinalTemp=LinksNotFailedFinal;
    NodesFailedScenILabelTemp = NodesCanFailOrigLabel(nodeTable.FailureState); % Failed nodes
    Resilience=zeros(length(DamagedComponentsTable.NodeStartId),1);
    ResiTotWithoutNewInterLinks = calculate_resi_given_links(LinksNotFailedFinal);
    for iStep=1:nSteps
        iStep
        if DamagedComponentsImportanceSort(iStep,1)~=DamagedComponentsImportanceSort(iStep,2) % the case of damaged node
            LinksNotFailedFinalTemp=vertcat(LinksNotFailedFinalTemp(:,:),DamagedComponentsImportanceSort(iStep,1:2));
            ResiTotTemp=calculate_resi_given_links(LinksNotFailedFinalTemp);
            Resilience(iStep)=ResiTotTemp;

        else % the case of damaged link
            NodesFailedScenILabelTemp1 = NodesFailedScenILabelTemp(NodesFailedScenILabelTemp~=DamagedComponentsImportanceSort(iStep,1)); % Failed nodes
    %       NodesFailedScenILabelTemp1 should be empty by the end.
            LinksStartNodesFailedLogical = ismember(LinksNotFailedScenI(:,1),NodesFailedScenILabelTemp1);
            LinksNotFailedFinalTemp1 = LinksNotFailedScenI(~LinksStartNodesFailedLogical,:);
            
            LinksFailedEndNodesLogical = ismember(LinksNotFailedFinalTemp1(:,2),NodesFailedScenILabelTemp1);
            NodesFailedScenILabelTemp=NodesFailedScenILabelTemp1;
            
            LinksNotFailedFinalTemp2 = LinksNotFailedFinalTemp1(~LinksFailedEndNodesLogical,:);
            
            LinksNotFailedFinalTemp3=vertcat(LinksNotFailedFinalTemp2,LinksNotFailedFinalTemp);
            LinksNotFailedFinalTemp=unique_matrix_rows(LinksNotFailedFinalTemp3);
            
    %       LinksNotFailedFinalTable=array2table(LinksNotFailedFinal
    %       WaterPowerNodeOperational = unique(reshape(LinksNotFailedFinal,[],1));
            % GraphAllInitial = graph(LinksNotFailedFinalTemp(:,1),LinksNotFailedFinalTemp(:,2));
            % plot(GraphAllInitial)

            % check connectivity
            [ResiTotTemp,ResiWaterTemp,ResiPowerTemp] = calculate_resi_given_links(LinksNotFailedFinalTemp);
%             ResiWater(iStep,1)=ResiWaterTemp;
%             ResiPower(iStep,1)=ResiPowerTemp;
            Resilience(iStep)=ResiTotTemp;

        end

    end
    
    
    Resilience=[ResiTotWithoutNewInterLinks; Resilience];
    
    fontsize_label=14;
    fontsize_tick=12;
    lWidth=1.25;
    
    plot(1:nSteps+1,Resilience,':bo','LineWidth',lWidth)
    ylim([min(Resilience)-0.05 1.05])
    
    xlabel('Time step','FontSize',fontsize_label,'FontWeight','bold')
    ylabel('Resilience','FontSize',fontsize_label,'FontWeight','bold')
    
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', fontsize_tick)
    % set(gcf, 'Color', 'None') % make the background transparent
    grid on;
    grid minor

    box off;
    ax2 = axes('Position',get(gca,'Position'),...
               'XAxisLocation','top',...
               'YAxisLocation','right',...
               'Color','none',...
               'XColor','k','YColor','k');
    set(ax2,'YTick', []);
    set(ax2,'XTick', []);
    box on;
    
    
        % links first
        % LinksFailedScenI = LinksFailedLogic(iScen,:); % identify failed links
        LinksNotFailedScenI = LinksOrig(~linkTable.FailureState,:); % remaining links

        % remove any edge that has failed nodes incident to it.       
        % NodesFailedScenILogic = NodesFailedLogic(iScen,:);
        NodesFailedScenILabel = NodesCanFailOrigLabel(nodeTable.FailureState); % Failed nodes

        LinksStartNodesFailedLogical = ismember(LinksNotFailedScenI(:,1), NodesFailedScenILabel);
        LinksFailedStartNodes = LinksNotFailedScenI(~LinksStartNodesFailedLogical,:);
        LinksFailedEndNodesLogical = ismember(LinksFailedStartNodes(:,2), NodesFailedScenILabel);
        LinksNotFailedFinal = LinksFailedStartNodes(~LinksFailedEndNodesLogical,:);
%         LinksNotFailedFinalTable=array2table(LinksNotFailedFinal
        WaterPowerNodeOperational = unique(reshape(LinksNotFailedFinal,[],1));
        % GraphAllInitial = graph(LinksNotFailedFinal(:,1),LinksNotFailedFinal(:,2));
        % plot(GraphAllInitial)

        % check connectivity
        [ResiTotTemp,ResiWaterTemp,ResiPowerTemp] = calculate_resi_given_links(LinksNotFailedFinal);
        ResiWater(iScen,1)=ResiWaterTemp;
        ResiPower(iScen,1)=ResiPowerTemp;
        ResiTot(iScen)=ResiTotTemp; % proportion of operational nodes of the i-th row
        disp(iScen)     
    end
    
%     ImportanceLinks=zeros(length(DamagedLinks(:,1)),1);
%     for iDamagedLink=1:length(DamagedLinks(:,1))
%         LinksNotFailedFinalTemp=LinksNotFailedFinal;
%         LinksNotFailedFinalTemp=vertcat(LinksNotFailedFinalTemp(:,:),DamagedLinks(iDamagedLink,:));
%         ResiTotTemp = calculate_resi_given_links(LinksNotFailedFinalTemp);
%         ImportanceLinks(iDamagedLink)=(ResiTotTemp-ResiTotWithNewInterLinks)/ResiTotWithNewInterLinks;
%     end
    
%     [ImportanceLinksDescend,ImportanceLinksRowIndex]=sort(ImportanceLinks,'descend');
%     DamagedLinksImportance=[DamagedLinks(ImportanceLinksRowIndex,:) ImportanceLinksDescend]; 
    
%     ImportanceNodes=zeros(length(DamagedNodesLabel),1);
    
    % 
%     for iDamagedNodes=1:length(DamagedNodesLabel)
%         LinksNotFailedFinalTemp=LinksNotFailedFinal;
%         % append those link whose start and end nodes have been restored
%                 
%         NodesFailedScenILabelTemp = NodesCanFailOrigLabel(nodeTable.FailureState & (NodesCanFailOrigLabel~=DamagedNodesLabel(iDamagedNodes))'); % Failed nodes
% 
%         LinksStartNodesFailedLogical = ismember(LinksNotFailedScenI(:,1), NodesFailedScenILabelTemp);
%         LinksFailedStartNodes = LinksNotFailedScenI(~LinksStartNodesFailedLogical,:);
%         LinksFailedEndNodesLogical = ismember(LinksFailedStartNodes(:,2), NodesFailedScenILabelTemp );
%         LinksNotFailedFinalTemp = LinksFailedStartNodes(~LinksFailedEndNodesLogical,:);
%         
%         ResiTotTemp = calculate_resi_given_links(LinksNotFailedFinalTemp);
%         ImportanceNodes(iDamagedNodes)=(ResiTotTemp-ResiTotWithNewInterLinks)/ResiTotWithNewInterLinks;
%     end
%     
    
%     [ImportanceNodesDescend,ImportanceNodesRowIndex]=sort(ImportanceNodes,'descend');
%     DamagedNodesImportance=[DamagedNodesLabel(ImportanceNodesRowIndex) ImportanceNodesDescend];  
%     
%     for iStep=1:nSteps
%         % links first
%         % LinksFailedScenI = LinksFailedLogic(iScen,:); % identify failed links
%         LinksNotFailedScenI = LinksOrig(~linkTable.FailureState,:); % remaining links
% 
%         % remove any edge that has failed nodes incident to it       
%         % NodesFailedScenILogic = NodesFailedLogic(iScen,:);
%         NodesFailedScenILabel = NodesCanFailOrigLabel(nodeTable.FailureState); % Failed nodes
% 
%         LinksStartNodesFailedLogical = ismember(LinksNotFailedScenI(:,1), NodesFailedScenILabel);
%         LinksFailedStartNodes = LinksNotFailedScenI(~LinksStartNodesFailedLogical,:);
%         LinksFailedEndNodesLogical = ismember(LinksFailedStartNodes(:,2), NodesFailedScenILabel);
%         LinksNotFailedFinal = LinksFailedStartNodes(~LinksFailedEndNodesLogical,:);
% %         LinksNotFailedFinalTable=array2table(LinksNotFailedFinal
%         WaterPowerNodeOperational = unique(reshape(LinksNotFailedFinal,[],1));
%         % GraphAllInitial = graph(LinksNotFailedFinal(:,1),LinksNotFailedFinal(:,2));
%         % plot(GraphAllInitial)
% 
%         % check connectivity
%         [ResiTotTemp,ResiWaterTemp,ResiPowerTemp] = calculate_resi_given_links(LinksNotFailedFinal);
%         ResiWater(iScen,1)=ResiWaterTemp;
%         ResiPower(iScen,1)=ResiPowerTemp;
%         ResiTot(iScen)=ResiTotTemp; % proportion of operational nodes of the i-th row
%         disp(iScen)     
%     end

    % !!!write a function. Inputs: Link table. Return: resilience
    % of individual systems and total resilience.

    % mean_func_node_noInteraction = sum(functionality_xi(:,:),1)/num_simu;

    % mean_func_node_Interaction = sum(Functionality_xi(:,:),1)/numScen;
    % % save ('mean_func_node_Interaction2','mean_func_node_Interaction');
    % function_nodeType = zeros(size(demand_nodeIDs,1),2);
    % function_nodeType(:,2) = nodeTable.NodeID(demand_nodeIDs);
    % originally WaterNodeTypeID, but should use the initial numbering scheme WaterNodeID
    % instead of the numbering in Hazus represented by WaterNodeTypeID.
    ResiTotInitialTemp = mean(ResiTot);
    
    
end
    
% time_execution=round(etime(clock,time_start)/60,1);
% 
% clc;
% fprintf('The execution time for %d \n simulation runs is  %d mins.\n', [numSimu time_execution])
% 
% functinaliryRatio_ouput=[function_nodeType(:,2) mean(mean_func_node)'];
% xlswrite('functinaliryRatio.xlsx',functinaliryRatio_ouput);
% difference_prob_ScenarioSimulationsAndFailure=sum(double(NodesFailedLogic))/num_scenarios-probFailure_nodes(1,:);
% 

% %% figures
% % functionality by node
% figure('name','Mean Functionality by Node');
% hold on
% % load waterNodeID_adachi.mat
% load mean_func_node_Interaction2.mat
% load mean_func_node_noInteraction2.mat
% stem_without=stem(demand_nodeIDs,mean_func_node_noInteraction,'LineWidth',1);
% stem_with=stem(demand_nodeIDs,mean_func_node_Interaction,'LineWidth',1);

% diff2=mean_func_node_noInteraction-mean_func_node_Interaction;

% save 

% ax = gca;
% ax.YGrid = 'on';
% ax.GridLineStyle = '-';

% ylim([0 1]);
% xlabel('Node at pipe intersections','fontweight','bold','fontsize',14)
% ylabel('Mean functionality ratio','fontweight','bold','fontsize',14)
% xticks(demand_nodeIDs);
% xt = get(gca, 'XTick');
% set(gca, 'FontSize', 13)
% grid on;
% 
% box off;
% ax2 = axes('Position',get(gca,'Position'),...
%            'XAxisLocation','top',...
%            'YAxisLocation','right',...
%            'Color','none',...
%            'XColor','k','YColor','k');
% set(ax2,'YTick', []);
% set(ax2,'XTick', []);
% box on;
% legend({'Without interdependence','With interdependence'},'location','Best');
% 
% 
% % serviceability ratio
% 
% %     % histogram
% %     fig_service_histogram=figure('name','Histogram of the Mean Serviceability under Earthquake');
% %     hist(mean_service)
% %     xlabel('Mean Serviceability')
% %     ylabel('Frequncy')
% 
%     % pdf
%     fig_service_PDF=figure('name','PDF of the Mean Serviceability udner Earthquake');
%     [pdf_meanService,pdf_meanService_xi] = ksdensity(mean_service);
%     plot(pdf_meanService_xi,pdf_meanService,'LineWidth',3)
%     xlabel('Mean serviceability ratio','fontweight','bold','fontsize',12)
%     ylabel('Probability density','fontweight','bold','fontsize',12)
%     xt = get(gca, 'XTick');
%     set(gca, 'FontSize', 11)
%     % set(gcf, 'Color', 'None') % make the background transparent
%     grid on;
% 
%     box off;
%     ax2 = axes('Position',get(gca,'Position'),...
%                'XAxisLocation','top',...
%                'YAxisLocation','right',...
%                'Color','none',...
%                'XColor','k','YColor','k');
%     set(ax2,'YTick', []);
%     set(ax2,'XTick', []);
%     box on;
% %     saveas(gcf,fig_service_PDF,'epsc')
%     
% 
% 
% %----------------------------------------------------------
% % calculate PGA at nodes
% 
% clear;
% clc;
% 
% waterPowerNode=xlsread('waterPowerNode.xlsx');
% 
% epicenter = [35.3 -90.3];
% M_w = 7.7;
% 
% dist_km= distdim(distance(epicenter(1), epicenter(2), waterPowerNode(:,4), waterPowerNode(:,5)), 'deg', 'km');
%     
%     term1 = 3.79;
%     term2 = (0.298*(M_w - 6));
%     term3 = (0.0536*(M_w - 6)^2);
%     term4 = log10(dist_km);
%     term5 = 0.00135*dist_km;
%     
% PGA = 10.^(term1 + term2 - term3 -term4 - term5); % pga
% 
% 
% waterPowerNode(:,6) = dist_km;
% waterPowerNode(:,7) = PGA;
% 
% waterPowerNodeTable=array2table(waterPowerNode,'VariableNames',{'NodeID','NodeType','NodeIDLabel','Lat','Long','distToCenter','PGA'});
% 
% 
% % number of nodes, only the first 15 nodes can fail, though.
% % num_of_nodes = waterPowerNodeTable.NodeIDLabel < 3;
% num_of_nodes = waterPowerNodeTable.NodeIDLabel ~= NaN; % All can fail
% NodeIdCanFail = waterPowerNodeTable.NodeID(num_of_nodes);
% 
% gravityAcceleration=9.80665;
% PGA_in_g=PGA/gravityAcceleration/100;  
% 
% % lambda_PfByPGA = [log(1.5) log(1.5) log(3.0) log(0.47) log(0.7) log(0.9)];   % Table 1, T. Adachi, B.R. Ellingwood / Reliability Engineering and System Safety 93 (2008), page 84     
% % zeta_PfByPGA = [0.8 0.6 0.1 0.4 0.4 0.4]; 
% lambdaPGA = [log(1.5) log(2.0) log(5.0) log(1.2) log(1.3) log(1.4)];  
% % assumed to generate more reasonable; pumping stations have power backup 
% zetaPGA = [0.8 0.6 0.2 0.4 0.4 0.4]; 
% 
% waterPowerNodeTable.NodeType(waterPowerNodeTable.NodeType>100)=waterPowerNodeTable.NodeType(waterPowerNodeTable.NodeType>100)-97;
% 
% for iNodeCanFail = 1:size(NodeIdCanFail,1)
%     iNodeType = waterPowerNodeTable.NodeType(iNodeCanFail);
%     NodePf=normcdf(PGA_in_g(iNodeCanFail),lambdaPGA(iNodeType),zetaPGA(iNodeType)); 
%     waterPowerNodeTable.NodePf(iNodeCanFail) = NodePf;    
% end
% 
% % save('waterPowerNodeTable', 'waterPowerNodeTable');
% % clear;
% % clc
% % load waterPowerNodeTable.mat
% 
% 
% 
% %----------------------------------------------------
% % calculate pf of water and power lines
% 
% % Calculate the distance between interdependent nodes
% % clear;
% % clc;
% % set directory 
% % read data
% % header:LinkID	NodeStartID	NodeEndID	TotLength	StartNodeX	StartNodeY	EndNodeX	EndNodeY
% waterPowerLink=xlsread('waterPowerLink.xlsx');
% 
% waterPowerLink(148:164,4)=distdim(distance(waterPowerLink(148:164,5),waterPowerLink(148:164,6),...
%     waterPowerLink(148:164,7),waterPowerLink(148:164,8)), 'deg', 'feet');
% 
% WaterPowerLinkTable=array2table(waterPowerLink,'VariableNames',{'LinkID','NodeStartID','NodeEndID',...
%     'TotLength','StartNodeLat','StartNodeLong','EndNodeLat','EndNodeLong'});
% 
% 
% % calcualte PGV
% epicenter = [35.3 -90.3];
% M_w = 7.7;
% 
% StartNodeDist_km = distdim(distance(epicenter(1), epicenter(2), WaterPowerLinkTable.StartNodeLat, WaterPowerLinkTable.StartNodeLong), 'deg', 'km');
% EndNodeDist_km = distdim(distance(epicenter(1), epicenter(2), WaterPowerLinkTable.EndNodeLat, WaterPowerLinkTable.EndNodeLong), 'deg', 'km');
%      
% term1 = 2.04;
% term2 = (0.422*(M_w - 6));
% term3 = (0.0373*(M_w - 6)^2);
% 
% StartNodeTerm4 = log10(StartNodeDist_km);
% EndNodeTerm4 = log10(EndNodeDist_km);
%     
% WaterPowerLinkTable.StartNodePGV_ins= 10.^(term1 + term2 - term3 -StartNodeTerm4)./2.54; 
% WaterPowerLinkTable.EndNodePGV_ins = 10.^(term1 + term2 - term3 -StartNodeTerm4)./2.54;     
% 
% 
% % water links
% a_waterlink=0.004;%0.00187;
% K_waterLink=0.5;
% 
% waterLinkID=1:71;
% RRWaterLink=a_waterlink.*K_waterLink.*(WaterPowerLinkTable.StartNodePGV_ins(waterLinkID)+...
%     WaterPowerLinkTable.EndNodePGV_ins(waterLinkID))./2;
% 
% WaterPowerLinkTable.LinkPf(waterLinkID) = 1-exp(-RRWaterLink.*WaterPowerLinkTable.TotLength(waterLinkID)./1000);
% 
% % power lines
% 
% a_powerlink=a_waterlink/2;
% K_powerLink=K_waterLink;
% 
% powerLinkID=72:147;
% RR_powerLink=a_powerlink.*K_powerLink.*(WaterPowerLinkTable.StartNodePGV_ins(powerLinkID)+...
%     WaterPowerLinkTable.EndNodePGV_ins(powerLinkID))./2;
% 
% WaterPowerLinkTable.LinkPf(powerLinkID) = 1-exp(-RR_powerLink.*WaterPowerLinkTable.TotLength(powerLinkID)./1000);
% 
% % interdependent links: water nodes demand electric power supply. 148,155
% % a and K is the same as in power lines
% 
% 
% 
% 
% 
% %-------------------------------------
% NodeTable = waterPowerNodeTable;
% LinkTable = WaterPowerLinkTable;
% 
% num_demandNodes=34;
% 
% % create base case u-v pairs: the ID of starting and end nodes respectively
% s_base = LinkTable.NodeStartID(:);
% t_base = LinkTable.NodeEndID(:);
% 
% numSimu= 100; 
% 
% nodes_can_Fail_IDs = NodeTable.NodeType ~= 3;
% nodes_can_Fail = NodeTable.NodeID(nodes_can_Fail_IDs);
% 
% % demand_nodeLogic = nodeTable.NodeType == 3 | nodeTable.NodeType == 5| nodeTable.NodeType == 6; % end-user water nodes and 23 and 12kv substations
% demand_nodeLogic = NodeTable.NodeType == 3; 
% demand_nodeIDs = NodeTable.NodeID(demand_nodeLogic); %16-49
% 
% num_nodes = size(nodes_can_Fail,1);
% num_links = size(LinkTable,1);
% 
% Ist_power=[0 1 2 4 6 8];
% mean_func_node_Interaction=zeros(length(Ist_power),num_demandNodes);
% 
% for i=1:length(Ist_power)
% 
%     interLinkPowerSupply=148:155; % supply water to power stations
%     a_interLinkPowerSupply=a_powerlink.*1;
%     RR_interLinkPowerSupply=a_interLinkPowerSupply.*K_powerLink.*(LinkTable.StartNodePGV_ins(interLinkPowerSupply)+...
%         LinkTable.EndNodePGV_ins(interLinkPowerSupply))./2;
% 
%     LinkTable.LinkPf(interLinkPowerSupply) = 1-exp(-RR_interLinkPowerSupply.*LinkTable.TotLength(interLinkPowerSupply)./1000);
%     LinkTable.LinkPf(148:155)=1;  % disable water to power link
%     % interdependent links: power nodes demand water supply. 156:164
% 
%     interLinkWaterSupply=156:164; % supply power to pumping stations
%     a_interLinkWaterSupply=a_waterlink.*Ist_power(i);
%     RR_interLinkWaterSupply=a_interLinkWaterSupply.*K_powerLink.*(LinkTable.StartNodePGV_ins(interLinkWaterSupply)+...
%         LinkTable.EndNodePGV_ins(interLinkWaterSupply))./2;
% 
%     LinkTable.LinkPf(interLinkWaterSupply) = 1-exp(-RR_interLinkWaterSupply.*LinkTable.TotLength(interLinkWaterSupply)./1000);
% 
%       
%     
%     U_nodes = rand(numSimu, num_nodes);
%     U_links = rand(numSimu, num_links);% generate num_scenarios*num_links vectors of random #
% 
%     probFailure_nodes = repmat(NodeTable.NodePf(nodes_can_Fail_IDs)',numSimu,1);
%     probFailure_links = repmat(LinkTable.LinkPf',numSimu,1);
%    
%     
%     % if U < pf_fail, then the node/link is disabled 
%     NodesFailedLogic = U_nodes <=probFailure_nodes;  % acceptance-rejection sampling to generate pf*num_scenarios
%     LinksFailedLogic = U_links <=probFailure_links;
%        
% 
%     supply_Logic = NodeTable.NodeType == 4 | NodeTable.NodeType == 5 | NodeTable.NodeType == 6; % one directional input dependence. water pumps depend on eletricity
%     supply_nodeIDs = NodeTable.NodeID(supply_Logic);
%     Functionality_xi = zeros(numSimu, size(demand_nodeIDs,1));
% 
% 
% NodeTable = waterPowerNodeTable;
% LinkTable = WaterPowerLinkTable;
% 
% 
% % num_simu= 5000; 
% % create base case u-v pairs: the ID of starting and end nodes respectively
% s_base = LinkTable.NodeStartID(:);
% t_base = LinkTable.NodeEndID(:);
% 
% nodes_can_Fail_IDs = NodeTable.NodeType ~= 3;
% nodes_can_Fail = NodeTable.NodeID(nodes_can_Fail_IDs);
% 
% % demand_nodeLogic = nodeTable.NodeType == 3 | nodeTable.NodeType == 5| nodeTable.NodeType == 6; % end-user water nodes and 23 and 12kv substations
% demand_nodeLogic = NodeTable.NodeType == 3; 
% demand_nodeIDs = NodeTable.NodeID(demand_nodeLogic); %16-49
% 
% num_nodes = size(nodes_can_Fail,1);
% num_links = size(LinkTable,1);
% 
% mean_func_node_noInteraction=zeros(length(Ist_power),num_demandNodes);
% 
% for i=1:length(Ist_power)
%     interLinkPowerSupply=148:155;
%     a_interLinkPowerSupply=a_powerlink*1;
%     RR_interLinkPowerSupply=a_interLinkPowerSupply.*K_powerLink.*(LinkTable.StartNodePGV_ins(interLinkPowerSupply)+...
%         LinkTable.EndNodePGV_ins(interLinkPowerSupply))./2;
% 
%     LinkTable.LinkPf(interLinkPowerSupply) = 1-exp(-RR_interLinkPowerSupply.*LinkTable.TotLength(interLinkPowerSupply)./1000);
% 
%     % interdependent links: power nodes demand water supply. 156:164
% 
%     interLinkWaterSupply=156:164;
%     a_interLinkWaterSupply=a_waterlink*Ist_power(i);
%     RR_interLinkWaterSupply=a_interLinkWaterSupply.*K_powerLink.*(LinkTable.StartNodePGV_ins(interLinkWaterSupply)+...
%         LinkTable.EndNodePGV_ins(interLinkWaterSupply))./2;
% 
%     LinkTable.LinkPf(interLinkWaterSupply) = 1-exp(-RR_interLinkWaterSupply.*LinkTable.TotLength(interLinkWaterSupply)./1000);
% 
%     %-----------------------
%     % a node in water distribution network is connected to an electrc power
%     % node after removing the dependence link: electric depends on water links
% 
%     % eletricNodeDemandWaterID=201:208;
%          for iScenario = 1:numSimu   % control+i and shift+tab
%              % links first
%              LinksOrig = [s_base t_base];        % original links
%              LinksFailedScenI = LinksFailedLogic(iScenario,:); % identify failed links
%              LinksNotFailedScenI = LinksOrig(~LinksFailedScenI',:); % remove links that failed in a simulation
% 
%              % then the nodes
%              NodesCanFailOrigLabel = nodes_can_Fail';       
%              NodesFailedScenILogic = NodesFailedLogic(iScenario,:);
%              NodesNotFailedScenILabel = NodesCanFailOrigLabel(NodesFailedScenILogic);
% 
%              % remove any instance of that node from edge list.(Remove failed nodes???)
%              NodesFailedStartLogical = ismember(LinksNotFailedScenI(:,1) , NodesNotFailedScenILabel);
%              LinksFailedStartNodes = LinksNotFailedScenI(~NodesFailedStartLogical,:);
%              LinksFailedEndNodesLogical = ismember(LinksFailedStartNodes(:,2) , NodesNotFailedScenILabel);
%              LinksNotFailedFinal = LinksFailedStartNodes(~LinksFailedEndNodesLogical,:);
% 
%              % check connectivity
%              G_i = graph();
%              G_i = addnode(G_i, max(nodeTable.NodeLabel));
%              G_i = addedge(G_i, LinksNotFailedFinal(:,1), LinksNotFailedFinal(:,2));
% 
%              DistToNodesInLoop = distances(G_i);
% 
%              num_isNotInf = sum(double(~isinf(DistToNodesInLoop)),1);
%              
%              Functionality_xi(iScenario,:) = double(~num_isNotInf == 0);   % takes on 1 or 0. 1 is counted as functionable.
%              disp(iScenario)
%          end
% 
%     mean_func_node_noInteraction(i,:) = sum(Functionality_xi(:,:),1)/numSimu;
% 
% end
% 
% diff=mean_func_node_noInteraction-mean_func_node_Interaction;
% 
% % set(0,'defaultaxeslinestyleorder',{'+','o','*','.','x','s','d','^'});
% 
% % % mrk=['ko-' 'gs-' 'b*' 'c+' 'm^' 'rv'];
% % for i=1:(length(Ist_power)-1)
% %     plot(demand_nodeIDs, diff(i,:))
% % end
% % legend('No Interaction','Ist0','2*Ist0','4*Ist0','6*Ist0','location','Best');
% 
% hold on
% % plot(demand_nodeIDs, diff(1,:),'ko-')
% plot(demand_nodeIDs, diff(2,:),'bs-')
% plot(demand_nodeIDs, diff(3,:),'r*--')
% plot(demand_nodeIDs, diff(4,:),'kv-')
% plot(demand_nodeIDs, diff(5,:),'bp-')
% plot(demand_nodeIDs, diff(6,:),'r+-')
% legend('Ist0','2*Ist0','4*Ist0','6*Ist0','location','Best');
% 
% 
% ax = gca;
% ax.YGrid = 'on';
% ax.GridLineStyle = '-';
% 
% ylim([0 0.6]);
% xlabel('Intermediate water delivery node','fontweight','bold','fontsize',12)
% ylabel('Difference between mean functionality ratios','fontweight','bold','fontsize',12)
% xticks(demand_nodeIDs);
% xt = get(gca, 'XTick');
% set(gca, 'FontSize', 13)
% 
% box off;
% ax2 = axes('Position',get(gca,'Position'),...
%            'XAxisLocation','top',...
%            'YAxisLocation','right',...
%            'Color','none',...
%            'XColor','k','YColor','k');
% set(ax2,'YTick', []);
% set(ax2,'XTick', []);
% box on;
% % legend([stem_without,stem_with],{'Without physical interdependence','With physical interdependence'},'location','Best');
% 
% 
% %_____________________
% % functionality ratio vs link length
% hold off
% 
% plot(LinkTable.TotLength(demand_nodeIDs),'r+-')
% 
% 
% % difference VS dependence strength
% 
% 
% 
% % hold on
% % stem(demand_nodeIDs, diff(1,:))
% % stem(demand_nodeIDs, diff(2,:))
% 
% 
% % serviceability ratio
% 
% %     % histogram
% %     fig_service_histogram=figure('name','Histogram of the Mean Serviceability under Earthquake');
% %     hist(mean_service)
% %     xlabel('Mean Serviceability')
% %     ylabel('Frequncy')
% 
%     % pdf
%     fig_service_PDF=figure('name','PDF of the Mean Serviceability udner Earthquake');
%     [pdf_meanService,pdf_meanService_xi] = ksdensity(mean_service);
%     plot(pdf_meanService_xi,pdf_meanService,'LineWidth',3)
%     xlabel('Mean serviceability ratio','fontweight','bold','fontsize',12)
%     ylabel('Probability density','fontweight','bold','fontsize',12)
%     xt = get(gca, 'XTick');
%     set(gca, 'FontSize', 11)
%     % set(gcf, 'Color', 'None') % make the background transparent
%     grid on;
% 
%     box off;
%     ax2 = axes('Position',get(gca,'Position'),...
%                'XAxisLocation','top',...
%                'YAxisLocation','right',...
%                'Color','none',...
%                'XColor','k','YColor','k');
%     set(ax2,'YTick', []);
%     set(ax2,'XTick', []);
%     box on;
% %     saveas(gcf,fig_service_PDF,'epsc')
% 
% 
% % Untitled 5
% 
%     i=2
%     interLinkPowerSupply=148:155; % supply water to power stations
%     a_interLinkPowerSupply=a_powerlink.*1;
%     RR_interLinkPowerSupply=a_interLinkPowerSupply.*K_powerLink.*(LinkTable.StartNodePGV_ins(interLinkPowerSupply)+...
%         LinkTable.EndNodePGV_ins(interLinkPowerSupply))./2;
% 
%     LinkTable.LinkPf(interLinkPowerSupply) = 1-exp(-RR_interLinkPowerSupply.*LinkTable.TotLength(interLinkPowerSupply)./1000);
%     LinkTable.LinkPf(148:155)=1;  % disable water to power link
%     % interdependent links: power nodes demand water supply. 156:164
% 
%     interLinkWaterSupply=156:164; % supply power to pumping stations
%     a_interLinkWaterSupply=a_waterlink.*Ist_power(i);
% %     RR_interLinkWaterSupply=a_interLinkWaterSupply.*K_powerLink.*(LinkTable.StartNodePGV_ins(interLinkWaterSupply)+...
% %         LinkTable.EndNodePGV_ins(interLinkWaterSupply))./2;
% 
% RR_interLinkWaterSupply=[0.01 0.02 0.04 0.06];
% hold on
% Ist=zeros(length(RR_interLinkWaterSupply),length(interLinkWaterSupply));
% for i=1:length(RR_interLinkWaterSupply)
%     Ist(i,:) = 1-exp(-RR_interLinkWaterSupply(i).*LinkTable.TotLength(interLinkWaterSupply)./1000);
%    
%     
% %     Ist=LinkTable.LinkPf(interLinkWaterSupply);
% 
% %     plot(sort(LinkTable.TotLength(interLinkWaterSupply)/1000),sort(LinkTable.LinkPf(interLinkWaterSupply)),'-d')
% end
%  
% % 
% hold on
% % plot(sort(LinkTable.TotLength(interLinkWaterSupply/1000),sort(LinkTable.LinkPf(interLinkWaterSupply)),'ko-')
% plot(sort(LinkTable.TotLength(interLinkWaterSupply)./1000), sort(Ist(1,:)),'ks-');
% plot(sort(LinkTable.TotLength(interLinkWaterSupply)./1000), sort(Ist(2,:)),'r*-');
% plot(sort(LinkTable.TotLength(interLinkWaterSupply)./1000), sort(Ist(3,:)),'bv-');
% plot(sort(LinkTable.TotLength(interLinkWaterSupply)./1000), sort(Ist(4,:)),'mv-');
% legend('I0','2*I0','4*I0','6*I0','location','Best');
% 
% 
% hold on
% % plot(sort(LinkTable.TotLength(interLinkWaterSupply/1000),sort(LinkTable.LinkPf(interLinkWaterSupply)),'ko-')
% plot(LinkTable.TotLength(interLinkWaterSupply)./1000, Ist(1,:),'ks-');
% plot(LinkTable.TotLength(interLinkWaterSupply)./1000, Ist(2,:),'r*-');
% plot(LinkTable.TotLength(interLinkWaterSupply)./1000, Ist(3,:),'bv-');
% plot(LinkTable.TotLength(interLinkWaterSupply)./1000, Ist(4,:),'mv-');
% legend('I0','2*I0','4*I0','6*I0','location','Best');
% 
% ax = gca;
% ax.YGrid = 'on';
% ax.GridLineStyle = '-';
% 
% ylim([0 0.6]);
% xlabel('Distance (1000 feet)','fontweight','bold','fontsize',12)
% ylabel('Interdependence strength','fontweight','bold','fontsize',12)
% 
% box off;
% ax2 = axes('Position',get(gca,'Position'),...
%            'XAxisLocation','top',...
%            'YAxisLocation','right',...
%            'Color','none',...
%            'XColor','k','YColor','k');
% set(ax2,'YTick', []);
% set(ax2,'XTick', []);
% box on;
% grid on
% legend([stem_without,stem_with],{'Without physical interdependence','With physical interdependence'},'location','Best');


