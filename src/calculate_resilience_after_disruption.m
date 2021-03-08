% function ResiTotInitialTemp=calculate_resilience_after_disruption(WaterPowerNodeTableTemp, WaterPowerLinkArray,nScen)
%%%%%%%%
%WaterPowerLinkTable=PowerAndWaterAllLinkSample(:,:,iSampleInterLink);
%%%%%%%%%

function [ResilTotStaticAll,ResilTotDynamicAll] = calculate_resilience_after_disruption(nScen,nStep,Mw)
    load WaterPowerNodeTableTemp
    load WaterPowerLinkArray
    %%%%%%%%
%     WaterPowerLinkArray = PowerAndWaterAllLinkSample(:,:,iSampleInterLink);
    %%%%%%%%%

    % node table: NodeLabel, NodeType, Lat, Long.
    % link table: link id, start and end node id, coordinates.
    epicenter = [35.3 -90.3];
%     Mw = 7.7;

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

    % % include population of each node
    % WaterPowerNodeTableTemp.Popul(:)=zeros(1,height(WaterPowerNodeTableTemp));
    % WaterPowerNodeTableTemp.Popul(WaterPowerNodeTableTemp.NodeLabel<100)=WaterNodeAttr(:,8);
    % WaterPowerNodeTableTemp.Popul(WaterPowerNodeTableTemp.NodeLabel>100)=PowerNodeAttr(:,11);


    % calcualte PGV for links

    %%%%%%%%%%%
    %  WaterPowerLinkArray=PowerAndWaterAllLinkSample(:,:,iSampleInterLink);
    %%%%%%%%%%%% 

    WaterPowerLinkTable=array2table(WaterPowerLinkArray,'VariableNames',...
        {'NodeStartId','NodeEndId','TotLength','StartNodeLat',...
        'StartNodeLong','EndNodeLat','EndNodeLong'});

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
    linkTable = WaterPowerLinkTable;

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

%     rng;
    U_nodes = rand(nScen, num_nodes);
%     rng;
    U_links = rand(nScen, num_links);% generate num_scenarios*num_links vectors of random #

    probFailure_nodes = repmat(nodeTable.NodePf(nodes_can_Fail_IDs)',nScen,1);
    probFailure_links = repmat(linkTable.LinkPf',nScen,1);

    % if U < pf_fail, then the node/link is disabled 
    NodesFailedLogic = U_nodes < probFailure_nodes;  % acceptance-rejection sampling to generate pf*num_scenarios
%     sum(NodesFailedLogic(1,:))
    LinksFailedLogic = U_links < probFailure_links;

    % Functionality_xi = zeros(numScen, max(nodeTable.NodeLabel));
    % nSteps=5;
    % TimeSteps= linspace(1,nSteps,nSteps);

    % ResiTot = zeros(nScen, nSteps);
    % ResiWater=zeros(nScen, nSteps);
    % ResiPower=zeros(nScen, nSteps);

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

%   nStep=80;
%   nStep=nStep;
    ResilTotDynamicAll=ones(nScen,nStep);
    ResilTotStaticAll=ones(nScen,nStep);

    % iScen = 1;
    for iScen = 1:nScen
%        fprintf('iteration of scenario: %d.\n', iScen)
      
        % initial resilience

        linkTable.FailureState=LinksFailedLogic(iScen,:)';
        nodeTable.FailureState=NodesFailedLogic(iScen,:)';

        % find failed power node supply to water
        PowerFailToSupplyWaterNode_Label=nodeTable.NodeLabel(nodeTable.NodeInterdependType...
                                        =="PowerSupplyWater" & nodeTable.FailureState==1);

        WaterFailToSupplyPowerNode_Label=nodeTable.NodeLabel(nodeTable.NodeInterdependType...
                                        =="WaterSupplyPower" & nodeTable.FailureState==1);

        %%%                            
        %%%%  WaterPowerNodeOperational
            % links first
        % LinksFailedScenI = LinksFailedLogic(iScen,:); % identify failed links
        LinksNotFailedScenIMayHaveFailedNodes = LinksOrig(~linkTable.FailureState,:); % remaining links

        % remove any edge that has failed nodes incident to it.       
        % NodesFailedScenILogic = NodesFailedLogic(iScen,:);
        NodesFailedScenILabel = NodesCanFailOrigLabel(nodeTable.FailureState); % Failed nodes

        LinksStartNodesFailedLogic = ismember(LinksNotFailedScenIMayHaveFailedNodes(:,1), NodesFailedScenILabel);
        LinksFailedStartNodes = LinksNotFailedScenIMayHaveFailedNodes(~LinksStartNodesFailedLogic,:);
        LinksFailedEndNodesLogic = ismember(LinksFailedStartNodes(:,2), NodesFailedScenILabel);
        LinksNotFailedFinalNoFailedNodes = LinksFailedStartNodes(~LinksFailedEndNodesLogic,:);
    %         LinksNotFailedFinalTable=array2table(LinksNotFailedFinal
        WaterPowerNodeOperational = unique(reshape(LinksNotFailedFinalNoFailedNodes,[],1));


        %%%%%                            

        % find the correpsonding water node that demands power

        WaterNodeNeedNewPowerNode_Label=linkTable.NodeEndId(ismember(linkTable.NodeStartId,...
                                        PowerFailToSupplyWaterNode_Label) & linkTable.NodeEndId<100);
        PowerNodeNeedNewWaterNode_Label=linkTable.NodeEndId(ismember(linkTable.NodeStartId,...
                                        WaterFailToSupplyPowerNode_Label) & linkTable.NodeEndId>100);

        PowerNodeOperational=nodeTable.NodeLabel(ismember(nodeTable.NodeLabel,...
                                        WaterPowerNodeOperational) & nodeTable.NodeType==6); %  end user power node type id=6
        PowerNodeOperationalRowLogic=ismember(nodeTable.NodeLabel,PowerNodeOperational);
        WaterNodeOperational=nodeTable.NodeLabel(ismember(nodeTable.NodeLabel,...
                                        WaterPowerNodeOperational) & nodeTable.NodeType==3); %  end user water node type id=3
        WaterNodeOperationalRowLogic=ismember(nodeTable.NodeLabel,WaterNodeOperational);

        NewPowerSupplyToWaterNodeLabel=zeros(length(WaterNodeNeedNewPowerNode_Label),1);
        NewWaterSupplyToPowerNodeLabel=zeros(length(PowerNodeNeedNewWaterNode_Label),1);
        % find the neareast alternative, operational (nodes in the final link table???) power supply node


        for iWaterNode=1:length(WaterNodeNeedNewPowerNode_Label)

            DistKmTemp= distdim(distance(nodeTable.Lat(iWaterNode), nodeTable.Long(iWaterNode), ...
                nodeTable.Lat(PowerNodeOperationalRowLogic),nodeTable.Long(PowerNodeOperationalRowLogic)), 'deg', 'km');
            if length(DistKmTemp)==1
                DistKmTempRevNorm=1;
            else
                % normalize distance
                normFactor=0.01;
                DistKmTempRev=1./DistKmTemp;
                range_DistKmTempRev = range(DistKmTempRev);
                DistKmTempRevNorm=(DistKmTempRev-min(DistKmTempRev)+normFactor*range_DistKmTempRev)./((1+normFactor)*range_DistKmTempRev);
            end
            ScorePowerNode=wtDist*DistKmTempRevNorm+wtPopul*nodeTable.PopulNorm(PowerNodeOperationalRowLogic)+...
                            wtSovi*nodeTable.Sovi(PowerNodeOperationalRowLogic);
            NewPowerSupplyToWaterNodeLabel(iWaterNode,1)=PowerNodeOperational(ScorePowerNode==max(ScorePowerNode));  

        end

        for iPowerNode=1:length(PowerNodeNeedNewWaterNode_Label)
            DistKmTemp= distdim(distance(nodeTable.Lat(iPowerNode), nodeTable.Long(iPowerNode), ...
                nodeTable.Lat(WaterNodeOperationalRowLogic),nodeTable.Long(WaterNodeOperationalRowLogic)), 'deg', 'km');
            if length(DistKmTemp)==1
                DistKmTempRevNorm=1;
            else
            % normalize distance
%             NormAdd=0.01;
                DistKmTempRev=1./DistKmTemp;
                range_DistKmTempRev = range(DistKmTempRev);
                DistKmTempRevNorm=(DistKmTempRev-min(DistKmTempRev)+normFactor*range_DistKmTempRev)./((1+normFactor)*range_DistKmTempRev);
            end

            ScoreWaterNode=wtDist*DistKmTempRevNorm+wtPopul*nodeTable.PopulNorm(WaterNodeOperationalRowLogic)+...
                            wtSovi*nodeTable.Sovi(WaterNodeOperationalRowLogic);
            NewWaterSupplyToPowerNodeLabel(iPowerNode,1)=WaterNodeOperational(ScoreWaterNode==max(ScoreWaterNode));  

        end

        InterLinksNew = vertcat([NewWaterSupplyToPowerNodeLabel,PowerNodeNeedNewWaterNode_Label],...
                       [NewPowerSupplyToWaterNodeLabel,WaterNodeNeedNewPowerNode_Label]);
        InterLinksNewTable=array2table(InterLinksNew, 'VariableNames',{'NodeStartId','NodeEndId'});
        InterLinksNewTable.FailureState=ones(length(InterLinksNewTable.NodeStartId),1);
        % New link table after building new interdependency links
        clear LinkTableWithNewInterLinks
        LinkTableWithNewInterLinks.NodeStartId = vertcat(linkTable.NodeStartId,InterLinksNewTable.NodeStartId);
        LinkTableWithNewInterLinks.NodeEndId = vertcat(linkTable.NodeEndId,InterLinksNewTable.NodeEndId);                 
        LinkTableWithNewInterLinks.FailureState = vertcat(linkTable.FailureState,InterLinksNewTable.FailureState); 

        LinkTableWithNewInterLinks=struct2table(LinkTableWithNewInterLinks);

    %     nSteps=sum(linkTable.FailureState)+sum(nodeTable.FailureState);
    %     TimeSteps = linspace(1,nSteps,nSteps);


        % component importance, will be used as the activity priority 

        % operational links, including new interlinks
        clear LinksNotFailedFinalWithNewInterLinks
        LinksNotFailedFinalWithNewInterLinks.NodeStartId = vertcat(LinksNotFailedFinalNoFailedNodes(:,1),...
                                                           InterLinksNewTable.NodeStartId);
        LinksNotFailedFinalWithNewInterLinks.NodeEndId = vertcat(LinksNotFailedFinalNoFailedNodes(:,2),...
                                                            InterLinksNewTable.NodeEndId);                 
        LinksNotFailedFinalWithNewInterLinks=struct2table(LinksNotFailedFinalWithNewInterLinks);

        % calculate the resilience after disruption with new inter links
        [ResilTotTemp,ResiWaterTemp,ResiPowerTemp] = calculate_resi_given_links(LinksNotFailedFinalWithNewInterLinks);
        ResiWaterWithNewInterLinks=ResiWaterTemp;
        ResiPowerWithNewInterLinks=ResiPowerTemp;
        ResiTotWithNewInterLinks=ResilTotTemp;


        % rank the importance of damaged links first

        % identify damaged links: failure state==1
        DamagedLinks.NodeStartId = vertcat(linkTable.NodeStartId(linkTable.FailureState==1),InterLinksNewTable.NodeStartId);
        DamagedLinks.NodeEndId = vertcat(linkTable.NodeEndId(linkTable.FailureState==1),InterLinksNewTable.NodeEndId);        
        % damagedLinks=struct2array(DamagedLinks);
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
        nDamagedComponents=length(DamagedComponentsTable.NodeStartId);
        ImportanceComponents=zeros(nDamagedComponents,1);
        ResilTotWithoutNewInterLinks = calculate_resi_given_links(LinksNotFailedFinalNoFailedNodes);

        for iDamagedComponent=1:length(DamagedComponentsTable.NodeStartId)
    %         iDamagedComponent
            LinksNotFailedFinalTemp1=LinksNotFailedFinalNoFailedNodes;
            if DamagedComponentsTable.ComponentType(iDamagedComponent)=="Link"
                LinksNotFailedFinalTemp1=vertcat(LinksNotFailedFinalTemp1(:,:),DamagedComponentsArray(iDamagedComponent,1:2));
                ResilTotTemp=calculate_resi_given_links(LinksNotFailedFinalTemp1);
                ImportanceComponents(iDamagedComponent)=(ResilTotTemp-ResilTotWithoutNewInterLinks)/ResilTotWithoutNewInterLinks;
            else
                NodesFailedScenILabelTemp1 = NodesCanFailOrigLabel(nodeTable.FailureState & (NodesCanFailOrigLabel~=...
                DamagedNodesArray(iDamagedComponent-height(DamagedLinksTable)))'); % remove the current node from table of failed nodes

                LinksStartNodesFailedLogic = ismember(LinksNotFailedScenIMayHaveFailedNodes(:,1), NodesFailedScenILabelTemp1);          
                LinksFailedEndNodesLogic = ismember(LinksNotFailedScenIMayHaveFailedNodes(:,2),NodesFailedScenILabelTemp1);
                LinksWithFailedNodesLogic=LinksStartNodesFailedLogic+LinksFailedEndNodesLogic;
                LinksNotFailedFinalTemp1 = LinksNotFailedScenIMayHaveFailedNodes(~LinksWithFailedNodesLogic,:);

                ResilTotTemp=calculate_resi_given_links(LinksNotFailedFinalTemp1);
                ImportanceComponents(iDamagedComponent)=(ResilTotTemp-ResilTotWithoutNewInterLinks)/ResilTotWithoutNewInterLinks;
        % length(LinksNotFailedFinalTemp(:,1))
        % Perhaps in calculating the importance, use resilience without new
        % interlinks.
            end

        end

        [ImportanceComponentsDescend,ImportanceComponentsRowIndex]=sort(ImportanceComponents,'descend');
        DamagedComponentsImportanceSort=[DamagedComponentsArray(ImportanceComponentsRowIndex,1:2) ImportanceComponentsDescend]; 
        DamagedComponentsImportanceSort;


        % restore the components and obtain the resilience curve
        nSteps=length(DamagedComponentsImportanceSort(:,1));
        LinksNotFailedFinalTemp = LinksNotFailedFinalNoFailedNodes;
        NodesFailedScenILabelTemp = NodesCanFailOrigLabel(nodeTable.FailureState); % Failed nodes
        Resil=zeros(length(DamagedComponentsTable.NodeStartId),1);
        ResilWater=zeros(length(DamagedComponentsTable.NodeStartId),1);
        ResilPower=zeros(length(DamagedComponentsTable.NodeStartId),1);
        for iStep=1:nSteps
%             fprintf('iteration of time steps: %d.\n',iStep)

            if DamagedComponentsImportanceSort(iStep,1)~=DamagedComponentsImportanceSort(iStep,2) % the case of damaged links
                LinksNotFailedFinalTemp=vertcat(LinksNotFailedFinalTemp(:,:),DamagedComponentsImportanceSort(iStep,1:2));
                [ResilTotTemp,ResiWaterTemp,ResiPowerTemp] = calculate_resi_given_links(LinksNotFailedFinalTemp);
                ResilWater(iStep)=ResiWaterTemp;
                ResilPower(iStep)=ResiPowerTemp;
                Resil(iStep)=ResilTotTemp;

            else % the case of damaged nodes
                NodesFailedScenILabelTemp = NodesFailedScenILabelTemp(NodesFailedScenILabelTemp~=DamagedComponentsImportanceSort(iStep,1)); % Failed nodes
            %   NodesFailedScenILabelTemp1 should be empty by the end.
                LinksStartNodesFailedLogic = ismember(LinksNotFailedScenIMayHaveFailedNodes(:,1),NodesFailedScenILabelTemp);
                LinksFailedEndNodesLogic = ismember(LinksNotFailedScenIMayHaveFailedNodes(:,2),NodesFailedScenILabelTemp);
                LinksWithFailedNodesLogic=LinksStartNodesFailedLogic+LinksFailedEndNodesLogic;
                LinksNotFailedFinalTemp1 = LinksNotFailedScenIMayHaveFailedNodes(~LinksWithFailedNodesLogic,:); 
                LinksNotFailedFinalTemp=union(LinksNotFailedFinalTemp1,LinksNotFailedFinalTemp,'rows');

        %       LinksNotFailedFinalTable=array2table(LinksNotFailedFinal
        %       WaterPowerNodeOperational = unique(reshape(LinksNotFailedFinal,[],1));
                % GraphAllInitial = graph(LinksNotFailedFinalTemp(:,1),LinksNotFailedFinalTemp(:,2));
                % plot(GraphAllInitial)

                % check connectivity
                [ResilTotTemp,ResiWaterTemp,ResiPowerTemp] = calculate_resi_given_links(LinksNotFailedFinalTemp);
                ResilWater(iStep)=ResiWaterTemp;
                ResilPower(iStep)=ResiPowerTemp;
                Resil(iStep)=ResilTotTemp;
            end
        end

        % dynamic componenet importance ranking
        ResilTotDynamic=zeros(nDamagedComponents,1);
        ResilWaterDynamic=zeros(nDamagedComponents,1);
        ResilPowerDynamic=zeros(nDamagedComponents,1);

        ImportanceRemainCompsDynamic=zeros(nDamagedComponents,6); % Initialization
        RemainDamagedCompsArray=DamagedComponentsArray;
        global LinksNotFailedFinalTemp ResilTotWithoutNewInterLinks NodesFailedScenILabelTemp LinksNotFailedScenIMayHaveFailedNodes
        LinksNotFailedFinalTemp=LinksNotFailedFinalNoFailedNodes; % initialize the links not fail table that changes each iter
        NodesFailedScenILabelTemp = NodesCanFailOrigLabel(nodeTable.FailureState); % Failed nodes
        for iComp=1:nDamagedComponents % 33 
    %         %   sort the components
    %         iComp     
    %     %   print('damaged components')
    %         RemainDamagedCompsArray
    %     %   print('link not failed final temp')
    %         LinksNotFailedFinalTemp
    %         length(LinksNotFailedFinalTemp(:,1)) 
            ImportanceRemainCompsSort=rank_components_importantce_dynamic(RemainDamagedCompsArray);
            ResilTotDynamic(iComp)=ImportanceRemainCompsSort(1,3);
            ResilWaterDynamic(iComp)=ImportanceRemainCompsSort(1,5);
            ResilPowerDynamic(iComp)=ImportanceRemainCompsSort(1,6);
            %   ResilienceDynamical(iComp)
            ImportanceRemainCompsDynamic(iComp,:)=ImportanceRemainCompsSort(1,:);        
            %   update remaining damaged comp list
            RemainDamagedCompsArray=ImportanceRemainCompsSort(2:end,1:2);
            %   update links not failed 
            if  ImportanceRemainCompsSort(1,1)~=ImportanceRemainCompsSort(1,2)
                LinksNotFailedFinalTemp=vertcat(LinksNotFailedFinalTemp(:,:),ImportanceRemainCompsSort(1,1:2));
            else % update failed nodes first
                % remove the current restored node from failed nodes list
                NodesFailedScenILabelTemp1 = NodesFailedScenILabelTemp(NodesFailedScenILabelTemp~=ImportanceRemainCompsSort(1,1)); 
                LinksStartNodesFailedLogic = ismember(LinksNotFailedScenIMayHaveFailedNodes(:,1),NodesFailedScenILabelTemp1);         
                LinksFailedEndNodesLogic = ismember(LinksNotFailedScenIMayHaveFailedNodes(:,2),NodesFailedScenILabelTemp1);
                LinksWithFailedNodesLogic=LinksStartNodesFailedLogic+LinksFailedEndNodesLogic;
                LinksNotFailedFinalTemp1 = LinksNotFailedScenIMayHaveFailedNodes(~LinksWithFailedNodesLogic,:);
                % add links that are restored due to restored nodes
                LinksNotFailedFinalTemp=union(LinksNotFailedFinalTemp1,LinksNotFailedFinalTemp,'rows');
                NodesFailedScenILabelTemp=NodesFailedScenILabelTemp1;
            end
        end

        ResilTotStatic=[ResilTotWithoutNewInterLinks; Resil];
        ResilTotDynamic=[ResilTotWithoutNewInterLinks; ResilTotDynamic];

        ResilTotStaticAll(iScen,1:length(ResilTotStatic))=ResilTotStatic';

        ResilTotDynamicAll(iScen,1:length(ResilTotDynamic))=ResilTotDynamic';
        
            
    end 

    if sum(ResilTotDynamicAll(:,length(ResilTotDynamic))==0)==nScen
        disp(sum(ResilTotDynamicAll(:,length(ResilTotDynamic))==0))
        msg = 'The last column is equal to a zero vec.\n';
        error(msg)
    end
end
       
