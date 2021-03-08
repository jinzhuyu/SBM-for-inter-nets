clc;
clear all;

% main part of code. Execute first.

%---------------------------------------------
% calculate probability of physical interlinks

% load data
% water node id:1-49; power: 157-117
WaterPowerNode = csvread('waterPowerNodeNew.csv',1,0);
WaterPowerNodeTable = array2table(WaterPowerNode,'VariableNames',...
    {'NodeLabel','NodeType','Lat','Long'});

% include population
PowerNodeAttr = xlsread('powerNodeAttribute.xlsx');

CencusPopul = xlsread('cencusPopulation.xlsx');
PowerNodeAttr(:,11) = zeros(length(PowerNodeAttr(:,7)),1);

for ii = 1:length(PowerNodeAttr(:,7))
    row_same = find(CencusPopul(:,1)==round(PowerNodeAttr(ii,7),0));
    PowerNodeAttr(ii,11) = CencusPopul(row_same,2);
end

WaterNodeAttr = xlsread('WaterNodeAttributes.xlsx');

WaterNodeAttr(:,8) = zeros(length(WaterNodeAttr(:,5)),1);

for ii = 1:length(WaterNodeAttr(:,5))
    row_same = find(CencusPopul(:,1)==round(WaterNodeAttr(ii,5),0));
    WaterNodeAttr(ii,8) = CencusPopul(row_same,2);
end

% include population of each node
WaterPowerNodeTable.Popul(:) = zeros(1,height(WaterPowerNodeTable));
WaterPowerNodeTable.Popul(WaterPowerNodeTable.NodeLabel<100) = WaterNodeAttr(:,8);
WaterPowerNodeTable.Popul(WaterPowerNodeTable.NodeLabel>100) = PowerNodeAttr(:,11);

WaterPowerNodeTable.Sovi(:) = zeros(1,height(WaterPowerNodeTable));
WaterPowerNodeTable.Sovi(WaterPowerNodeTable.NodeLabel<100) = WaterNodeAttr(:,7);
WaterPowerNodeTable.Sovi(WaterPowerNodeTable.NodeLabel>100) = PowerNodeAttr(:,9);

% normalize population
normLowBound= 0.01;  % can be other values, such as 0.005. Meaning of a constant should be given.
WaterPowerNodeTable.PopulNorm(:) = zeros(1,height(WaterPowerNodeTable));
WaterPowerNodeTable.PopulNorm = (WaterPowerNodeTable.Popul-min(WaterPowerNodeTable.Popul)+normLowBound*range(WaterPowerNodeTable.Popul))./...
                                ((1+normLowBound)*range(WaterPowerNodeTable.Popul));
                            
% figure()
% hist(WaterPowerNodeTable.PopulNorm)
% probability of interlink from power to water between type

% intermidiate water delivery node (type id:3) supply to power gate stations (type id:101)
% select NodeType==3
WaterNodeInterDelivery = WaterPowerNodeTable.NodeLabel(WaterPowerNodeTable.NodeType==3);
% select NodeType==101
PowerNodeGateStation = WaterPowerNodeTable.NodeLabel(WaterPowerNodeTable.NodeType==101);

nPowerNodeGateStation = length(PowerNodeGateStation);
nWaterNodeInterDelivery = length(WaterNodeInterDelivery);
PrWater2Power = zeros(nPowerNodeGateStation,nWaterNodeInterDelivery);
nRow = nPowerNodeGateStation*nWaterNodeInterDelivery;
nCol = 8;
PrWater2PowerTable = array2table(zeros(nRow,nCol),'VariableNames',...
                        {'NodeStartId','NodeEndId','TotLength','StartNodeLat',...
                        'StartNodeLong','EndNodeLat','EndNodeLong','Pr'});

wtDist = 0.5;
wtPopul = 0.25;
wtSovi = 1-wtDist-wtPopul;
distFromWaterNodesKm = zeros(nPowerNodeGateStation,nWaterNodeInterDelivery);

for iGateStation = 1:nPowerNodeGateStation
    GateStationId = PowerNodeGateStation(iGateStation);
    GateStationRowId = WaterPowerNodeTable.NodeLabel==GateStationId;
%     distFromWaterNodesKm = zeros(1,nWaterNodeInterDelivery);
    for jInterDelivery = 1:nWaterNodeInterDelivery
        InterDeliveryId = WaterNodeInterDelivery(jInterDelivery);
        InterDeliveryRowId = WaterPowerNodeTable.NodeLabel==InterDeliveryId;
        distFromWaterNodesKm(iGateStation,jInterDelivery) = ...
            distdim( distance( WaterPowerNodeTable.Lat(GateStationRowId),...
                               WaterPowerNodeTable.Long(GateStationRowId),...
                               WaterPowerNodeTable.Lat(InterDeliveryRowId),...
                               WaterPowerNodeTable.Long(InterDeliveryRowId) ),...
                               'deg','km');
    end
    
    % normalize distance
    normLowBound = 0.01;
    distRecipFromWaterKmTemp = 1./distFromWaterNodesKm(iGateStation,:);
    distRecipFromWaterKmTempNorm = truncated_normalize(distRecipFromWaterKmTemp,normLowBound);
%     figure()
%     hist(distKmTempRecipNorm)
    for jInterDelivery = 1:nWaterNodeInterDelivery
        InterDeliveryId = WaterNodeInterDelivery(jInterDelivery);
        InterDeliveryRowId = WaterPowerNodeTable.NodeLabel==InterDeliveryId;
        % use find
           
        PrWater2Power(iGateStation,jInterDelivery) = wtDist*distRecipFromWaterKmTempNorm(jInterDelivery)+...
                                                    wtPopul*WaterPowerNodeTable.PopulNorm(InterDeliveryRowId)+...
                                                    wtSovi*WaterPowerNodeTable.Sovi(InterDeliveryRowId);
                                                
        IndexTableRow = (iGateStation-1)*nWaterNodeInterDelivery+jInterDelivery;
        PrWater2PowerTable.NodeStartId(IndexTableRow) = InterDeliveryId;
        PrWater2PowerTable.NodeEndId(IndexTableRow) = GateStationId;
        PrWater2PowerTable.TotLength(IndexTableRow) = distFromWaterNodesKm(jInterDelivery);
        PrWater2PowerTable.StartNodeLat(IndexTableRow) = WaterPowerNodeTable.Lat(InterDeliveryRowId);
        PrWater2PowerTable.StartNodeLong(IndexTableRow) = WaterPowerNodeTable.Long(InterDeliveryRowId);
        PrWater2PowerTable.EndNodeLat(IndexTableRow) = WaterPowerNodeTable.Lat(GateStationRowId);
        PrWater2PowerTable.EndNodeLong(IndexTableRow) = WaterPowerNodeTable.Long(GateStationRowId);
    end 
    % normalization
    PrWater2Power(iGateStation,:) = PrWater2Power(iGateStation,:)./sum(PrWater2Power(iGateStation,:));
end
% reshape into a column vector
PrWater2PowerColVec = reshape(PrWater2Power',[],1);
% include Start node, end node, Start node X (Lat), Start node Y (Long), end node X, end
% node Y
PrWater2PowerTable.Pr = PrWater2PowerColVec;


% intermidiate 12kv substations (type id:103) supply to pumping stations (type id:2)
% select NodeType==3
WaterNodePumpStation = WaterPowerNodeTable.NodeLabel(WaterPowerNodeTable.NodeType==2);
% select NodeType==101
% select 12 kv substation
PowerNode12Substation = WaterPowerNodeTable.NodeLabel(WaterPowerNodeTable.NodeType==103);

nPowerNode12Substation = length(PowerNode12Substation);
nWaterNodePumpStation = length(WaterNodePumpStation);
nRow = nWaterNodePumpStation*nPowerNode12Substation;
nCol = 8;
PrPower2WaterTable = array2table(zeros(nRow,nCol),'VariableNames',...
    {'NodeStartId','NodeEndId','TotLength','StartNodeLat',...
    'StartNodeLong','EndNodeLat','EndNodeLong','Pr'});

PrPower2Water = zeros(nWaterNodePumpStation,nPowerNode12Substation);
distFromPowerNodesKm = zeros(nWaterNodePumpStation,nPowerNode12Substation);
for iPumpStation = 1:nWaterNodePumpStation
    PumpStationId = WaterNodePumpStation(iPumpStation);
    PumpStationRowIndex = WaterPowerNodeTable.NodeLabel==PumpStationId;
    
    for j12Substation = 1:nPowerNode12Substation
        SubstationId = PowerNode12Substation(j12Substation);
        SubstationRowIndex = WaterPowerNodeTable.NodeLabel==SubstationId;
        distFromPowerNodesKm(iPumpStation,j12Substation) = ...
            distdim( distance( WaterPowerNodeTable.Lat(PumpStationRowIndex),...
                               WaterPowerNodeTable.Long(PumpStationRowIndex),...
                               WaterPowerNodeTable.Lat(SubstationRowIndex),...
                               WaterPowerNodeTable.Long(SubstationRowIndex) ),...
                               'deg', 'km' );
    end
    
    % normalize distance
    normLowBound = 0.01;
    distRecipFromPowerKmTemp = 1./distFromPowerNodesKm(iPumpStation,:);
    distRecipFromPowerKmTempNorm = truncated_normalize(distRecipFromPowerKmTemp,normLowBound);
    
    for j12Substation = 1:nPowerNode12Substation
        SubstationId = PowerNode12Substation(j12Substation);
        SubstationRowIndex = WaterPowerNodeTable.NodeLabel==SubstationId;
        % probability of interdependency links

        PrPower2Water(iPumpStation,j12Substation) = wtDist*distRecipFromPowerKmTempNorm(j12Substation)+...
                                                    wtPopul*WaterPowerNodeTable.PopulNorm(SubstationRowIndex)+...
                                                    wtSovi*WaterPowerNodeTable.Sovi(SubstationRowIndex);
  
        IndexTableRow = (iPumpStation-1)*nPowerNode12Substation+j12Substation;
        PrPower2WaterTable.NodeStartId(IndexTableRow) = SubstationId;
        PrPower2WaterTable.NodeEndId(IndexTableRow) = PumpStationId;
        PrPower2WaterTable.TotLength(IndexTableRow) = distFromWaterNodesKm(j12Substation);
        PrPower2WaterTable.StartNodeLat(IndexTableRow) = WaterPowerNodeTable.Lat(SubstationRowIndex);
        PrPower2WaterTable.StartNodeLong(IndexTableRow) = WaterPowerNodeTable.Long(SubstationRowIndex);
        PrPower2WaterTable.EndNodeLat(IndexTableRow) = WaterPowerNodeTable.Lat(PumpStationRowIndex);
        PrPower2WaterTable.EndNodeLong(IndexTableRow) = WaterPowerNodeTable.Long(PumpStationRowIndex);
    end 
    % normalize the probability
    PrPower2Water(iPumpStation,:) = PrPower2Water(iPumpStation,:)./sum(PrPower2Water(iPumpStation,:));
end
% reshape into a column vector
PrPower2WaterColVec = reshape(PrPower2Water',[],1);

PrPower2WaterTable.Pr = PrPower2WaterColVec;

PrWaterAndPowerTable = vertcat(PrWater2PowerTable,PrPower2WaterTable);
% water nodeIlabel match link Endpoints's label. NodeObjectId==NodeLabel.
% power nodeId not equal to objectId. So the NodeObjectiveId is switched to NodeLabel


%---------------------------------
% sampling interdependent networks

% Import link table. First without interlinks.
WaterPowerLink = csvread('waterPowerLinkNoInterlinks.csv',1,0);
WaterPowerLinkTable = array2table(WaterPowerLink,'VariableNames',...
    {'LinkId','NodeStartId','NodeEndId','TotLength','StartNodeLat',...
    'StartNodeLong','EndNodeLat','EndNodeLong'});

% sampling interlinks and merge with links of separate networks
% preprocess links of separate networks
WaterPowerLinkTableTemp = removevars(WaterPowerLinkTable,'LinkId');
WaterPowerLinkTableTemp.TotLength = WaterPowerLinkTableTemp.TotLength/3280.84; % length in km
WaterPowerLinkArray = table2array(WaterPowerLinkTableTemp);

nSampleInterlink = 30;
Water2PowerInterlinkSample = zeros(nPowerNodeGateStation,nCol-1,nSampleInterlink);
Power2WaterInterlinkSample = zeros(nWaterNodePumpStation,nCol-1,nSampleInterlink);
PowerAndWaterInterlinkSample = zeros(nWaterNodePumpStation+nPowerNodeGateStation,nCol-1,nSampleInterlink);
nPowerAndWaterAllLinkSample = nWaterNodePumpStation+nPowerNodeGateStation+...
                              length(WaterPowerLinkTable.NodeStartId);
PowerAndWaterAllLinkSample = zeros(nPowerAndWaterAllLinkSample,nCol-1,nSampleInterlink);
for iSampleInterlink = 1:nSampleInterlink
    %selcte a the probabilities for a Start node
    for jGateStation = 1:nPowerNodeGateStation
        RowIndex = ((jGateStation-1)*nWaterNodeInterDelivery+1):(jGateStation*nWaterNodeInterDelivery);
        PrRows = PrWater2PowerTable.Pr(RowIndex)';
        % sample one link
        PrRowsSample = gendist(PrRows,1,1); % pRows should be row vector. Returns the indice (RowIndexSample).
        RowsSampleIndex = (jGateStation-1)*nWaterNodeInterDelivery+PrRowsSample;             
        Water2PowerInterlinkSample(jGateStation,:,iSampleInterlink) = PrWater2PowerTable{RowsSampleIndex,1:(nCol-1)};
    end

    %selcte a the probabilities for a Start node
    for jPumpStation = 1:nWaterNodePumpStation
        RowIndex = ((jPumpStation-1)*nPowerNode12Substation+1):(jPumpStation*nPowerNode12Substation);
        PrRows = PrPower2WaterTable.Pr(RowIndex)';
        % sample one link
        PrRowsSample = gendist(PrRows,1,1); % pRows should be row vector. Returns the indice (RowIndexSample).
        RowsSampleIndex = (jPumpStation-1)*nPowerNode12Substation+PrRowsSample;             
        Power2WaterInterlinkSample(jPumpStation,:,iSampleInterlink) = PrPower2WaterTable{RowsSampleIndex,1:(nCol-1)};
    end
    % concatenate arrays of links and interlinks
    PowerAndWaterInterlinkSample(:,:,iSampleInterlink) = vertcat(Water2PowerInterlinkSample(:,:,iSampleInterlink),...
        Power2WaterInterlinkSample(:,:,iSampleInterlink));
    % combine interlinks and links of separate networks    
    PowerAndWaterAllLinkSample(:,:,iSampleInterlink) = vertcat(WaterPowerLinkArray,...
                                                               Water2PowerInterlinkSample(:,:,iSampleInterlink),...
                                                               Power2WaterInterlinkSample(:,:,iSampleInterlink));
%     disp(iSampleInterlink)   
end

% resilience after disruption for each sample of interdependency links.
ResiTotAfterDisruption = zeros(nSampleInterlink,1);
WaterPowerNodeTableTemp = WaterPowerNodeTable;
save WaterPowerNodeTableTemp;
nScen = 20;
nStep = 90;
global WaterPowerNodeTableTemp WaterPowerLinkArray

% Mw_vec = (3:8)+0.5;
Mw_vec = 5:8;
nMw = length(Mw_vec);

nSampleInterlink = nSampleInterlink;
% iSampleInterlink = 5;

ResilStatic_all = zeros(nScen,nStep,nSampleInterlink,nMw);
ResilDynamic_all = zeros(nScen,nStep,nSampleInterlink,nMw);

% % ResilStatic_all_1=ResilStatic_all;
% % ResilDynamic_all_1=ResilDynamic_all;
% % save ResilStatic_all_1 ResilDynamic_all_1
% 
% ResilStatic_all_2 = cat(3,ResilStatic_all_1,ResilStatic_al);
% ResilDynamic_all_2 = cat(3,ResilDynamic_all_1,ResilDynamic_all);
% 
% save ResilStatic_all_2 ResilDynamic_all_2

% iSampleInterlink = 1;
n_err = 0;
% iMw = 1;
for iMw = 1:nMw
    iSampleInterlink = 1;
    while iSampleInterlink <= nSampleInterlink
        try
            fprintf('Iteration of samples of inter links: %d.\n', iSampleInterlink)
            WaterPowerLinkArray = PowerAndWaterAllLinkSample(:,:,iSampleInterlink);
            save WaterPowerLinkArray
            [A,B] = calculate_resilience_after_disruption(nScen,nStep,Mw_vec(iMw));
            ResilStatic_all(:,:,iSampleInterlink,iMw)= A;
            ResilDynamic_all(:,:,iSampleInterlink,iMw)= B;
            iSampleInterlink = iSampleInterlink+1;
        catch Err
            fprintf('Error: %s\n', Err.message);
    %         iSampleInterlink = iSampleInterlink;
            n_err = n_err + 1;
            if n_err > max(nSampleInterlink/10,3)
                disp(n_err)
                keyboard % stop the loop when excessive erros occur and 
                         % enter debugger mode to take care of the errors.

            else
                continue;
            end
        end
    end
end

save ResilDynamic_all ResilStatic_all
% generate graph from the link table
% for iSampleInterlink=1:nSampleInterlink
%     % generate graph from the link table
% %     sAllBase =PowerAndWaterAllLinkSample(:,1,iSample);
% %     tAllBase =PowerAndWaterAllLinkSample(:,2,iSample);
% % %   [uSBase, ~, uTBase] = unique([sBase tBase]);
% % %   NetworkInitial = graph(uTBase(1:end/2), uTBase(end/2+1:end), [], cellstr(num2str(uSBase.')));
% %     GraphAllInitial = graph(sAllBase,tAllBase);
% %     plot(GraphAllInitial)
%  
% %     figure('Name','water-power networks')
% %     plot(GraphAllInitial)
% %     AjacMat=adjacency(GraphAllInitial);
%     WaterPowerLinkArray=PowerAndWaterAllLinkSample(:,:,iSampleInterlink);
%     ResiTotAfterDisruption(iSampleInterlink,1)=calculate_resilience_after_disruption(WaterPowerNodeTableTemp, WaterPowerLinkArray,nScen);
% end

% % Estimate initial network flow by betweenness
% 
% % modify supply nodes and demand nodes
% 
%   % water networks
%     WaterNodeRowIndex=PowerAndWaterAllLinkSample(:,1,iSample)<100 & PowerAndWaterAllLinkSample(:,2,iSample)<100;
%     sWaterBase =PowerAndWaterAllLinkSample(WaterNodeRowIndex,1,iSample);
%     tWaterBase =PowerAndWaterAllLinkSample(WaterNodeRowIndex,2,iSample);
%     GraphWaterInitial = graph(sWaterBase,tWaterBase);
%  
% %     figure('Name','Water network')
% %     plot(GraphWaterInitial)
% %     xlim([-3.5 3.5])
% %     ylim([-3.5 3.5])
%     
%     WaterSupplyNode=WaterPowerNodeTable.NodeLabel(WaterPowerNodeTable.NodeType==1);
%     WaterDemandNode=WaterPowerNodeTable.NodeLabel(WaterPowerNodeTable.NodeType==3 | WaterPowerNodeTable.NodeType==2);
% 
%     % Note: have to find all shortest paths that traverse a given node.
%     % Concatenate shortest paths 
%     % Initialize the zero-th shortest path
%     kPaths=25; % set it to a relatively large number
%     WaterMaxDist=30;
% 
%     WaterNodeId=WaterPowerNodeTable.NodeLabel(WaterPowerNodeTable.NodeType<100);
%     nWaterNode=length(WaterNodeId);
% %     WaterShortPaths=zeros(1,WaterMaxDist);
%     nWaterSupplyNode=length(WaterSupplyNode);
%     nWaterDemandNode=length(WaterDemandNode);
%     WaterBCOD=zeros(nWaterSupplyNode,nWaterDemandNode,nWaterNode);
%     for iSupplyNode=1:nWaterSupplyNode
%         for jDemandNode=1:nWaterDemandNode
%         % ShortPathTemp=shortestpath(GraphWaterInitial,WaterSupplyNode(iSupplyNode),WaterDemandNode(jDemandNode));
%         [ShortPathsTemp,NumShortestPaths]=findallshortestpaths(GraphWaterInitial,kPaths,WaterSupplyNode(iSupplyNode),WaterDemandNode(jDemandNode));
% %         WaterShortPaths=padconcatenation(WaterShortPaths,ShortPathsTemp,1);
%             for kWaterNode=1:nWaterNode
% %                 fprintf('i:%d',iSupplyNode)
% %                 fprintf('j:%d', jDemandNode)
% %                 fprintf('k:%d',kPowerNode)
%                 WaterBCOD(iSupplyNode,jDemandNode,kWaterNode)=sum(sum(ShortPathsTemp==WaterNodeId(kWaterNode)))/NumShortestPaths;
%             
%             end
%         end
%         disp(iSupplyNode)
%     end
% 
%     % Count the frequancy of occurence for a given node that lies on a shortest path
% %     MaxLength=max(max(sWaterBase,tWaterBase));
% %     WaterFreqNodeInShortPath=ones(MaxLength,1);
% %     for iNode=1:MaxLength
% %         WaterFreqNodeInShortPath(iNode)=nnz(WaterShortPaths==iNode)';
% %     end
%     WaterBCODFinal=zeros(nWaterNode,1);
%     for kWaterNode=1:nWaterNode
%         WaterBCODFinal(kWaterNode)=sum(sum(WaterBCOD(:,:,kWaterNode)));
%     end
%     
%     
%  % Power flow data
%     % Estimate initial network flow by betweenness
%     % Power networks
%     PowerNodeRowIndex=PowerAndWaterAllLinkSample(:,1,iSample)>100 & PowerAndWaterAllLinkSample(:,2,iSample)>100;
%     sPowerBase =PowerAndWaterAllLinkSample(PowerNodeRowIndex,1,iSample);
%     tPowerBase =PowerAndWaterAllLinkSample(PowerNodeRowIndex,2,iSample);
%     GraphPowerInitial = graph(sPowerBase,tPowerBase);
%    
% %     figure('Name','Power network')
% %     plot(GraphPowerInitial)
% %     xlim([0 8])
% %     ylim([0 8])
%     
%     PowerSupplyNode=WaterPowerNodeTable.NodeLabel(WaterPowerNodeTable.NodeType==101);
%     PowerDemandNode=WaterPowerNodeTable.NodeLabel(WaterPowerNodeTable.NodeType==102 | WaterPowerNodeTable.NodeType==103);
% 
%     % note: have to find all shortest paths that traverse a given node.
%     % concatenate shortest paths
%     kPaths=15;
%     % Initialize the zero-th shortest path
% %     PowerMaxDist=25;
% %     PowerShortPaths=zeros(1,PowerMaxDist);
%     nPowerSupplyNode=length(PowerSupplyNode);
%     nPowerDemandNode=length(PowerDemandNode);
%     PowerNodeId=WaterPowerNodeTable.NodeLabel(WaterPowerNodeTable.NodeType>100);
%     nPowerNode=length(PowerNodeId);
% %     PowerNumShortestPaths=zeros(nPowerSupplyNode,nPowerDemandNode,nPowerNode);
% %     PowerNodeRowIndex=WaterPowerNodeTable.NodeLabel==PowerNodeId;
%     PowerBCOD=zeros(nPowerSupplyNode,nPowerDemandNode,nPowerNode);
%     for iSupplyNode=1:nPowerSupplyNode
% %         fprintf('i:%d',iSupplyNode)
%         for jDemandNode=1:nPowerDemandNode
%             
%             % ShortPathTemp=shortestpath(GraphPowerInitial,PowerSupplyNode(iSupplyNode),PowerDemandNode(jDemandNode));
%             [ShortPathsTemp,NumShortestPaths]=findallshortestpaths(GraphPowerInitial,kPaths,PowerSupplyNode(iSupplyNode),PowerDemandNode(jDemandNode));
% %             PowerShortPaths=padconcatenation(PowerShortPaths,ShortPathsTemp,1);
% %             PowerNumShortestPaths(iSupplyNode,jDemandNode)=NumShortestPaths;
%             for kPowerNode=1:nPowerNode
% %                 fprintf('i:%d',iSupplyNode)
% %                 fprintf('j:%d', jDemandNode)
% %                 fprintf('k:%d',kPowerNode)
%                 PowerBCOD(iSupplyNode,jDemandNode,kPowerNode)=sum(sum(ShortPathsTemp==PowerNodeId(kPowerNode)))/NumShortestPaths;
%             
%             end
%             
%         end
%    
%     end
% 
%     PowerBCODFinal=zeros(nPowerNode,1);
%     for kPowerNode=1:nPowerNode
%         PowerBCODFinal(kPowerNode)=sum(sum(PowerBCOD(:,:,kPowerNode)));
%     end    
%     
%     
%     % count the frequancy of occurence for a given node that lies on a shortest path
%     MaxLength=max(max(sPowerBase,tPowerBase));
%     PowerFreqNodeInShortPath=ones(MaxLength,1);
%     for iNode=1:MaxLength
%         PowerFreqNodeInShortPath(iNode)=nnz(PowerShortPaths==iNode)';
%     end

 
%   Intial resilience= portional of operational nodes
 
% for each samples of interdependent network
    % purtube the networks. Compute the failure probability of components
    
    
    % find direct damage nodes
    % find nodes inoperational due to cascading failure
    % reconstruct the graph
    % compute resilience
    
    % restore in 8 time steps
    % find those who are totally damaged and remove
    % restore based on capacity. restore 1/8 components
  
% Useful functions for finding multiple shortest paths

% %https://www.mathworks.com/matlabcentral/fileexchange/36086-modified-dij
% %sktra-s-algorithm-to-return-all-paths-that-tie-for-shortest
% 
% %ref.: https://www.mathworks.com/help/bioinfo/ref/graphshortestpath.html
% s=[1 2 3 4];
% t=[2 3 4 1];
% W =ones(1,4);
% UG=sparse(s,t,W);
% UG=tril(G+G');
% % view(biograph(UG,[],'ShowArrows','off','ShowWeights','on'))
% [dist,path,pred]=graphshortestpath(UG,3,1,'Directed',false)
% 
% dist 
% path
% %----------------
% % in real networks, multiple shortest paths are very unlikely
% % % use the function dijkstraties
% % weightMatrix=1e3*ones(length(s),length(t));
% % matriz_costo=weightMatrix;
% % [sp, spcost] = dijkstraties(matriz_costo, sIndex, dIndex);
