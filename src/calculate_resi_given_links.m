function [ResiTot,ResiWater,ResiPower] = calculate_resi_given_links(LinksTable)

% input: table of operational links: start node label, end node label 
% output: 1. total resilience. 2. Resilience of water and power networks
% respectively.
global WaterDemandNodeLabel;
global PowerDemandNodeLabel;
global nWaterDemandNode;
global nPowerDemandNode;

% if linksTable is a table variable, convert to array
if isa(LinksTable, 'table')==1
    LinksTable=table2array(LinksTable);
end

numNodeIdAll=max(max(LinksTable));

G_i = graph();
%  G_i = addnode(G_i, size(nodeTable.NodeID,1));
G_i = addnode(G_i, numNodeIdAll);
G_i = addedge(G_i, LinksTable(:,1), LinksTable(:,2));
%  figure(2); plot(G_i,'Layout','layered');

%  find edges (ID of incident Nodes) that completes a cycle. 
EdgesInLoops = dfsearch(G_i, 1, 'edgetodiscovered', 'Restart', true);

% select cycles that contain interlinks
EdgesInLoopsSelect=EdgesInLoops(floor(EdgesInLoops(:,1)./100)~=floor(EdgesInLoops(:,2)./100),:);
NodesInLoopId=EdgesInLoopsSelect(:);
NodesInLoopIdUnique=unique(NodesInLoopId);
DistToNodesInLoop = distances(G_i,NodesInLoopIdUnique);

num_isNotInf = sum(double(~isinf(DistToNodesInLoop)),1);
% count the number of not inf in each column sum(,1)
Functionality_xi = double(~num_isNotInf == 0); % takes on 1 or 0 with 1 representing functionable.

Weights=[0.5 0.5];
ResiWater=sum(Functionality_xi(WaterDemandNodeLabel))/nWaterDemandNode;
ResiPower=sum(Functionality_xi(PowerDemandNodeLabel))/nPowerDemandNode;
ResiTot=Weights*[ResiWater ResiPower]'; % proportion of operational nodes of the i-th row
   