% 2018 October. ICASP conference paper.
%----------------------------------------------------------
% calculate PGA at nodes

clear;
clc;

waterPowerNode=xlsread('waterPowerNode.xlsx');

epicenter = [35.3 -90.3];
M_w = 7.7;

dist_km= distdim(distance(epicenter(1), epicenter(2), waterPowerNode(:,4), waterPowerNode(:,5)), 'deg', 'km');
    
    term1 = 3.79;
    term2 = (0.298*(M_w - 6));
    term3 = (0.0536*(M_w - 6)^2);
    term4 = log10(dist_km);
    term5 = 0.00135*dist_km;
    
PGA = 10.^(term1 + term2 - term3 -term4 - term5); % pga


waterPowerNode(:,6) = dist_km;
waterPowerNode(:,7) = PGA;

waterPowerNodeTable=array2table(waterPowerNode,'VariableNames',{'NodeID','NodeType','NodeIDLabel','Lat','Long','distToCenter','PGA'});


% number of nodes, only the first 15 nodes can fail, though.
% num_of_nodes = waterPowerNodeTable.NodeIDLabel < 3;
num_of_nodes = waterPowerNodeTable.NodeIDLabel ~= NaN; % All can fail
NodeIDS_canFail = waterPowerNodeTable.NodeID(num_of_nodes);

gravityAcceleration=9.80665;
PGA_in_g=PGA/gravityAcceleration/100;  

% lambda_PfByPGA = [log(1.5) log(1.5) log(3.0) log(0.47) log(0.7) log(0.9)];   % Table 1, T. Adachi, B.R. Ellingwood / Reliability Engineering and System Safety 93 (2008), page 84     
% zeta_PfByPGA = [0.8 0.6 0.1 0.4 0.4 0.4]; 
lambda_PfByPGA = [log(1.5) log(2.0) log(2.5) log(1.2) log(1.3) log(1.4)];  
% assumed to generate more reasonable; pumping stations have power backup 
zeta_PfByPGA = [0.8 0.6 0.2 0.4 0.4 0.4]; 

waterPowerNodeTable.NodeType(waterPowerNodeTable.NodeType>100)=waterPowerNodeTable.NodeType(waterPowerNodeTable.NodeType>100)-97;

for index_i = 1:size(NodeIDS_canFail,1)
    i_type = waterPowerNodeTable.NodeType(index_i);
    NodePf=normcdf(PGA_in_g(index_i),lambda_PfByPGA(i_type),zeta_PfByPGA(i_type)); 
    waterPowerNodeTable.NodePf(index_i) = NodePf;    
end

% save('waterPowerNodeTable', 'waterPowerNodeTable');
% clear;
% clc
% load waterPowerNodeTable.mat



%----------------------------------------------------
% calculate pf of water and power lines

% Calculate the distance between interdependent nodes
clear;
clc;
% set directory 
% read data
% header:LinkID	NodeStartID	NodeEndID	TotLength	StartNodeX	StartNodeY	EndNodeX	EndNodeY
waterPowerLink=xlsread('waterPowerLink.xlsx');

waterPowerLink(148:164,4)=distdim(distance(waterPowerLink(148:164,5),waterPowerLink(148:164,6),...
    waterPowerLink(148:164,7),waterPowerLink(148:164,8)), 'deg', 'feet');

waterPowerLinkTable=array2table(waterPowerLink,'VariableNames',{'LinkID','NodeStartID','NodeEndID',...
    'TotLength','StartNodeLat','StartNodeLong','EndNodeLat','EndNodeLong'});


% calcualte PGV
epicenter = [35.3 -90.3];
M_w = 7.7;

StartNodeDist_km = distdim(distance(epicenter(1), epicenter(2), waterPowerLinkTable.StartNodeLat, waterPowerLinkTable.StartNodeLong), 'deg', 'km');
EndNodeDist_km = distdim(distance(epicenter(1), epicenter(2), waterPowerLinkTable.EndNodeLat, waterPowerLinkTable.EndNodeLong), 'deg', 'km');
     
term1 = 2.04;
term2 = (0.422*(M_w - 6));
term3 = (0.0373*(M_w - 6)^2);

StartNodeTerm4 = log10(StartNodeDist_km);
EndNodeTerm4 = log10(EndNodeDist_km);
    
waterPowerLinkTable.StartNodePGV_ins= 10.^(term1 + term2 - term3 -StartNodeTerm4)./2.54; 
waterPowerLinkTable.EndNodePGV_ins = 10.^(term1 + term2 - term3 -StartNodeTerm4)./2.54;     


% %create a nested structure for each pipe
% pipeSegments = struct('PipePoints', []);


% water links
a_waterlink=0.004;%0.00187;
K_waterLink=0.5;

waterLinkID=1:71;
RR_waterLink=a_waterlink.*K_waterLink.*(waterPowerLinkTable.StartNodePGV_ins(waterLinkID)+...
    waterPowerLinkTable.EndNodePGV_ins(waterLinkID))./2;

waterPowerLinkTable.LinkPf(waterLinkID) = 1-exp(-RR_waterLink.*waterPowerLinkTable.TotLength(waterLinkID)./1000);

% power lines

a_powerlink=a_waterlink/2;
K_powerLink=K_waterLink;

powerLinkID=72:147;
RR_powerLink=a_powerlink.*K_powerLink.*(waterPowerLinkTable.StartNodePGV_ins(powerLinkID)+...
    waterPowerLinkTable.EndNodePGV_ins(powerLinkID))./2;

waterPowerLinkTable.LinkPf(powerLinkID) = 1-exp(-RR_powerLink.*waterPowerLinkTable.TotLength(powerLinkID)./1000);

% interdependent links: water nodes demand electric power supply. 148,155
% a and K is the same as in power lines
interLinkPowerSupply=148:155;
a_interLinkPowerSupply=a_powerlink.*2; %%%%%%%%%%%%%%%%%%%%%!!!!!!!!!!!!!!!!!!!!!!
RR_interLinkPowerSupply=a_interLinkPowerSupply.*K_powerLink.*(waterPowerLinkTable.StartNodePGV_ins(interLinkPowerSupply)+...
    waterPowerLinkTable.EndNodePGV_ins(interLinkPowerSupply))./2;

waterPowerLinkTable.LinkPf(interLinkPowerSupply) = 1-exp(-RR_interLinkPowerSupply.*waterPowerLinkTable.TotLength(interLinkPowerSupply)./1000);

% interdependent links: power nodes demand water supply. 156:164

interLinkWaterSupply=156:164;
a_interLinkWaterSupply=a_waterlink.*1;
RR_interLinkWaterSupply=a_interLinkWaterSupply.*K_powerLink.*(waterPowerLinkTable.StartNodePGV_ins(interLinkWaterSupply)+...
    waterPowerLinkTable.EndNodePGV_ins(interLinkWaterSupply))./2;

waterPowerLinkTable.LinkPf(interLinkWaterSupply) = 1-exp(-RR_interLinkWaterSupply.*waterPowerLinkTable.TotLength(interLinkWaterSupply)./1000);

save('waterPowerLinkTable', 'waterPowerLinkTable');

% load waterPowerLinkTable.mat


%-------------------------------------------------
% Assess network functionality

clear all;

load waterPowerNodeTable.mat  
load waterPowerLinkTable.mat

nodeTable = waterPowerNodeTable;
linkTable = waterPowerLinkTable;


% create base case u-v pairs: the ID of starting and end nodes respectively
s_base = linkTable.NodeStartID(:);
t_base = linkTable.NodeEndID(:);

% generate grpah
G = graph();
G = addnode(G, size(nodeTable.NodeID,1));
G = addedge(G, s_base,t_base);
% plot(G)
% d = distances(G);

% generate u vectors, size of edges + links
% double loop Monte Carlo
num_simu= 2000; 
% scenarios of failure of nodes and links in the network. 
% tunr probilistic cases to deterministic scenarios: Pf is either 0 or 1. 

nodes_can_Fail_IDs = nodeTable.NodeType ~= 3;
nodes_can_Fail = nodeTable.NodeID(nodes_can_Fail_IDs);
% Note:The water distribution nodes are assumed to be undamaged because they are simply the branch points
% in the major water distribution system needed to distribute water to.(Takao Adachi,2009,pp.243)

% demand_nodeLogic = nodeTable.NodeType == 3 | nodeTable.NodeType == 5| nodeTable.NodeType == 6; % end-user water nodes and 23 and 12kv substations
demand_nodeLogic = nodeTable.NodeType == 3; 
demand_nodeIDs = nodeTable.NodeID(demand_nodeLogic); %16-49

% supply_Logic = nodeTable.NodeType == 2 | nodeTable.NodeType == 4; mutual
% independent
% supply_Logic = nodeTable.NodeType ==2; 

% interdependent. pumps depends on power. Supply nodes are power nodes.
supply_Logic = nodeTable.NodeType == 4 | nodeTable.NodeType == 5 | nodeTable.NodeType == 6; % one directional input dependence. water pumps depend on eletricity

supply_nodeIDs = nodeTable.NodeID(supply_Logic);
num_nodes = size(nodes_can_Fail,1);
num_links = size(linkTable,1);

func_node = zeros(num_simu,size(demand_nodeIDs,1));
service = zeros(num_simu,1);


% find nodes in a cycle

% calculate the distance to each of the nodes in the cycle

% if at least one distance is not Inf, then the node is in operation.

% iterations should be independent.Loop iterations must be consecutive, increasing integer values.
U_nodes = rand(num_simu, num_nodes);
U_links = rand(num_simu, num_links);% generate num_scenarios*num_links vectors of random #

probFailure_nodes = repmat(nodeTable.NodePf(nodes_can_Fail_IDs)',num_simu,1);
probFailure_links = repmat(linkTable.LinkPf',num_simu,1);

% if U < pf_fail, then the node/link is disabled 
f_nodes = U_nodes <probFailure_nodes;  % acceptance-rejection sampling to generate pf*num_scenarios
f_links = U_links <probFailure_links;

functionality_xi = zeros(num_simu, size(demand_nodeIDs,1));
serviceability_net = zeros(num_simu, 1);

%-----------------------
% a node in water distribution network is connected to an electrc power
% node after removing the dependence link: electric depends on water links

% eletricNodeDemandWaterID=201:208;
linkTable.LinkPf(148:155)=1;
 for i_scenario = 1:num_simu   % control+i and shift+tab
     % links first
     links_scen_orig = [s_base t_base];        % original links
     link_fail_scen_i = f_links(i_scenario,:); % identify failed links
     links_scen_i = links_scen_orig(~link_fail_scen_i',:); % remove links that failed in a simulation
     
     % then the nodes
     nodes_scen_i_orig = nodes_can_Fail';       
     node_fail_scen_i = f_nodes(i_scenario,:);
     nodes_scen_i = nodes_scen_i_orig(node_fail_scen_i);
     
     % remove any instance of that node from edge list.(Remove failed nodes???)
     remove_firstcolumn_logical = ismember(links_scen_i(:,1) , nodes_scen_i);
     links_scen_1 = links_scen_i(~remove_firstcolumn_logical,:);
     remove_secondcolumn_logical = ismember(links_scen_1(:,2) , nodes_scen_i);
     links_scen_final = links_scen_1(~remove_secondcolumn_logical,:);
     
     % check connectivity
     G_i = graph();
%    G_i = addnode(G_i, size(nodeTable.NodeID,1));
     G_i = addnode(G_i, max(nodeTable.NodeID));
     G_i = addedge(G_i, links_scen_final(:,1), links_scen_final(:,2));

     d_i = distances(G_i);
 
     % plot the graph
     %         h_i=plot(G_i);
     %         h_i.NodeColor='r';
     %         h_i.Marker='o';
    
     d_st=d_i(supply_nodeIDs, demand_nodeIDs); %pick out the distances between pairs of storage tanks and distribution nodes.
   
     num_isNotInf = sum(double(~isinf(d_st)),1);
     % count the number of not inf in each column sum(,1)
     % converts symbolic values to double precision
     functionality_xi(i_scenario,:) = double(~num_isNotInf == 0);   % takes on 1 or 0. 1 is counted as functionable.
     serviceability_net(i_scenario) = sum(functionality_xi(i_scenario,:),2)/size(demand_nodeIDs,1); % ratio of functioning node of the i-th row
     disp(i_scenario)
 end

%-----------------------
% there is a loop. For case of mutual depedence.
%  for i_scenario = 1:num_simu   % control+i and shift+tab
%      % links first
%      links_scen_orig = [s_base t_base];        % original links
%      link_fail_scen_i = f_links(i_scenario,:); % identify failed links
%      links_scen_i = links_scen_orig(~link_fail_scen_i',:); % remove links that failed in a simulation
%      
%      % then the nodes
%      nodes_scen_i_orig = nodes_can_Fail';       
%      node_fail_scen_i = f_nodes(i_scenario,:);
%      nodes_scen_i = nodes_scen_i_orig(node_fail_scen_i);
%      
%      % remove any instance of that node from edge list.(Remove failed nodes???)
%      remove_firstcolumn_logical = ismember(links_scen_i(:,1) , nodes_scen_i);
%      links_scen_1 = links_scen_i(~remove_firstcolumn_logical,:);
%      remove_secondcolumn_logical = ismember(links_scen_1(:,2) , nodes_scen_i);
%      links_scen_final = links_scen_1(~remove_secondcolumn_logical,:);
%      
%      % check connectivity
%      G_i = digraph();
% %      G_i = addnode(G_i, size(nodeTable.NodeID,1));
%      G_i = addnode(G_i, max(nodeTable.NodeID));
%      G_i = addedge(G_i, links_scen_final(:,1), links_scen_final(:,2));
% %      plot(G_i,'Layout','layered');
% 
%     % find edges that completes a cycle. % e contains the ID of Nodes in a
%     % cycle
%      [e] = dfsearch(G_i, 1, 'edgetodiscovered', 'Restart', true);
%      
%      ID_nodeInLoop=e(:);
%      ID_nodeInLoop_unique=unique(ID_nodeInLoop);
% %      d_i = distances(G_i);
%      d_i = distances(G_i,ID_nodeInLoop_unique);
%      %  Returns a matrix, d, where d(i,j) is the length of the shortest path between node i and node j.
%      %  If the graph is unweighted, all edge distances are taken to be 1.
%      
%      % plot the graph
%      %         h_i=plot(G_i);
%      %         h_i.NodeColor='r';
%      %         h_i.Marker='o';
%      
% %      d_st=d_i(supply_nodeIDs, demand_nodeIDs); %pick out the distances between pairs of storage tanks and distribution nodes.
%      d_st=d_i(:, demand_nodeIDs); %pick out the distances between pairs of storage tanks and distribution nodes.
%      
%      num_isNotInf = sum(double(~isinf(d_st)),1);
%      % count the number of not inf in each column sum(,1)
%      % converts symbolic values to double precision
%      functionality_xi(i_scenario,:) = double(~num_isNotInf == 0);   % takes on 1 or 0. 1 is counted as functionable.
%      serviceability_net(i_scenario) = sum(functionality_xi(i_scenario,:),2)/size(demand_nodeIDs,1); % ratio of functioning node of the i-th row
%      disp(i_scenario)
%  end
    
% mean_func_node_noInteraction = sum(functionality_xi(:,:),1)/num_simu;
% save ('mean_func_node_noInteraction','mean_func_node_noInteraction');


mean_func_node_Interaction = sum(functionality_xi(:,:),1)/num_simu;
save ('mean_func_node_Interaction2','mean_func_node_Interaction');
% function_nodeType = zeros(size(demand_nodeIDs,1),2);
% function_nodeType(:,2) = nodeTable.NodeID(demand_nodeIDs);
% originally WaterNodeTypeID, but should use the initial numbering scheme WaterNodeID
% instead of the numbering in Hazus represented by WaterNodeTypeID.
% mean_service = mean(serviceability_net);
    
% time_execution=round(etime(clock,time_start)/60,1);

% clc;
% fprintf('The execution time for %d \n simulation runs is  %d mins.\n', [num_simu time_execution])

% functinaliryRatio_ouput=[function_nodeType(:,2) mean(mean_func_node)'];
% xlswrite('functinaliryRatio.xlsx',functinaliryRatio_ouput);
% difference_prob_ScenarioSimulationsAndFailure=sum(double(f_nodes))/num_scenarios-probFailure_nodes(1,:);

%% Save the whole workspace without the figures to reduce the file size.
% save WaterNetworkEarthquakeHazard_Network_10012018.mat
% 
% pipeLinks.pf_constantPGV

% load CalculateNetworkServiceabilityUnderEarthquake.mat

%% figures
% functionality by node

figure('name','Mean Functionality by Node');

hold on
% load waterNodeID_adachi.mat
load mean_func_node_Interaction2.mat
load mean_func_node_noInteraction2.mat
stem_without=stem(demand_nodeIDs,mean_func_node_noInteraction,'LineWidth',1);
stem_with=stem(demand_nodeIDs,mean_func_node_Interaction,'LineWidth',1);


diff2=mean_func_node_noInteraction-mean_func_node_Interaction;

save 


ax = gca;
ax.YGrid = 'on';
ax.GridLineStyle = '-';

ylim([0 1]);
xlabel('Node at pipe intersections','fontweight','bold','fontsize',14)
ylabel('Mean functionality ratio','fontweight','bold','fontsize',14)
xticks(demand_nodeIDs);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 13)
grid on;

box off;
ax2 = axes('Position',get(gca,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k');
set(ax2,'YTick', []);
set(ax2,'XTick', []);
box on;
legend({'Without interdependence','With interdependence'},'location','Best');


% serviceability ratio

%     % histogram
%     fig_service_histogram=figure('name','Histogram of the Mean Serviceability under Earthquake');
%     hist(mean_service)
%     xlabel('Mean Serviceability')
%     ylabel('Frequncy')

    % pdf
    fig_service_PDF=figure('name','PDF of the Mean Serviceability udner Earthquake');
    [pdf_meanService,pdf_meanService_xi] = ksdensity(mean_service);
    plot(pdf_meanService_xi,pdf_meanService,'LineWidth',3)
    xlabel('Mean serviceability ratio','fontweight','bold','fontsize',12)
    ylabel('Probability density','fontweight','bold','fontsize',12)
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 11)
    % set(gcf, 'Color', 'None') % make the background transparent
    grid on;

    box off;
    ax2 = axes('Position',get(gca,'Position'),...
               'XAxisLocation','top',...
               'YAxisLocation','right',...
               'Color','none',...
               'XColor','k','YColor','k');
    set(ax2,'YTick', []);
    set(ax2,'XTick', []);
    box on;
%     saveas(gcf,fig_service_PDF,'epsc')
    


%----------------------------------------------------------
% calculate PGA at nodes

clear;
clc;

waterPowerNode=xlsread('waterPowerNode.xlsx');

epicenter = [35.3 -90.3];
M_w = 7.7;

dist_km= distdim(distance(epicenter(1), epicenter(2), waterPowerNode(:,4), waterPowerNode(:,5)), 'deg', 'km');
    
    term1 = 3.79;
    term2 = (0.298*(M_w - 6));
    term3 = (0.0536*(M_w - 6)^2);
    term4 = log10(dist_km);
    term5 = 0.00135*dist_km;
    
PGA = 10.^(term1 + term2 - term3 -term4 - term5); % pga


waterPowerNode(:,6) = dist_km;
waterPowerNode(:,7) = PGA;

waterPowerNodeTable=array2table(waterPowerNode,'VariableNames',{'NodeID','NodeType','NodeIDLabel','Lat','Long','distToCenter','PGA'});


% number of nodes, only the first 15 nodes can fail, though.
% num_of_nodes = waterPowerNodeTable.NodeIDLabel < 3;
num_of_nodes = waterPowerNodeTable.NodeIDLabel ~= NaN; % All can fail
NodeIDS_canFail = waterPowerNodeTable.NodeID(num_of_nodes);

gravityAcceleration=9.80665;
PGA_in_g=PGA/gravityAcceleration/100;  

% lambda_PfByPGA = [log(1.5) log(1.5) log(3.0) log(0.47) log(0.7) log(0.9)];   % Table 1, T. Adachi, B.R. Ellingwood / Reliability Engineering and System Safety 93 (2008), page 84     
% zeta_PfByPGA = [0.8 0.6 0.1 0.4 0.4 0.4]; 
lambda_PfByPGA = [log(1.5) log(2.0) log(5.0) log(1.2) log(1.3) log(1.4)];  
% assumed to generate more reasonable; pumping stations have power backup 
zeta_PfByPGA = [0.8 0.6 0.2 0.4 0.4 0.4]; 

waterPowerNodeTable.NodeType(waterPowerNodeTable.NodeType>100)=waterPowerNodeTable.NodeType(waterPowerNodeTable.NodeType>100)-97;

for index_i = 1:size(NodeIDS_canFail,1)
    i_type = waterPowerNodeTable.NodeType(index_i);
    NodePf=normcdf(PGA_in_g(index_i),lambda_PfByPGA(i_type),zeta_PfByPGA(i_type)); 
    waterPowerNodeTable.NodePf(index_i) = NodePf;    
end

% save('waterPowerNodeTable', 'waterPowerNodeTable');
% clear;
% clc
% load waterPowerNodeTable.mat



%----------------------------------------------------
% calculate pf of water and power lines

% Calculate the distance between interdependent nodes
% clear;
% clc;
% set directory 
% read data
% header:LinkID	NodeStartID	NodeEndID	TotLength	StartNodeX	StartNodeY	EndNodeX	EndNodeY
waterPowerLink=xlsread('waterPowerLink.xlsx');

waterPowerLink(148:164,4)=distdim(distance(waterPowerLink(148:164,5),waterPowerLink(148:164,6),...
    waterPowerLink(148:164,7),waterPowerLink(148:164,8)), 'deg', 'feet');

waterPowerLinkTable=array2table(waterPowerLink,'VariableNames',{'LinkID','NodeStartID','NodeEndID',...
    'TotLength','StartNodeLat','StartNodeLong','EndNodeLat','EndNodeLong'});


% calcualte PGV
epicenter = [35.3 -90.3];
M_w = 7.7;

StartNodeDist_km = distdim(distance(epicenter(1), epicenter(2), waterPowerLinkTable.StartNodeLat, waterPowerLinkTable.StartNodeLong), 'deg', 'km');
EndNodeDist_km = distdim(distance(epicenter(1), epicenter(2), waterPowerLinkTable.EndNodeLat, waterPowerLinkTable.EndNodeLong), 'deg', 'km');
     
term1 = 2.04;
term2 = (0.422*(M_w - 6));
term3 = (0.0373*(M_w - 6)^2);

StartNodeTerm4 = log10(StartNodeDist_km);
EndNodeTerm4 = log10(EndNodeDist_km);
    
waterPowerLinkTable.StartNodePGV_ins= 10.^(term1 + term2 - term3 -StartNodeTerm4)./2.54; 
waterPowerLinkTable.EndNodePGV_ins = 10.^(term1 + term2 - term3 -StartNodeTerm4)./2.54;     


% water links
a_waterlink=0.004;%0.00187;
K_waterLink=0.5;

waterLinkID=1:71;
RR_waterLink=a_waterlink.*K_waterLink.*(waterPowerLinkTable.StartNodePGV_ins(waterLinkID)+...
    waterPowerLinkTable.EndNodePGV_ins(waterLinkID))./2;

waterPowerLinkTable.LinkPf(waterLinkID) = 1-exp(-RR_waterLink.*waterPowerLinkTable.TotLength(waterLinkID)./1000);

% power lines

a_powerlink=a_waterlink/2;
K_powerLink=K_waterLink;

powerLinkID=72:147;
RR_powerLink=a_powerlink.*K_powerLink.*(waterPowerLinkTable.StartNodePGV_ins(powerLinkID)+...
    waterPowerLinkTable.EndNodePGV_ins(powerLinkID))./2;

waterPowerLinkTable.LinkPf(powerLinkID) = 1-exp(-RR_powerLink.*waterPowerLinkTable.TotLength(powerLinkID)./1000);

% interdependent links: water nodes demand electric power supply. 148,155
% a and K is the same as in power lines





%-------------------------------------
NodeTable = waterPowerNodeTable;
LinkTable = waterPowerLinkTable;

num_demandNodes=34;

% create base case u-v pairs: the ID of starting and end nodes respectively
s_base = LinkTable.NodeStartID(:);
t_base = LinkTable.NodeEndID(:);

num_simu= 10000; 

nodes_can_Fail_IDs = NodeTable.NodeType ~= 3;
nodes_can_Fail = NodeTable.NodeID(nodes_can_Fail_IDs);

% demand_nodeLogic = nodeTable.NodeType == 3 | nodeTable.NodeType == 5| nodeTable.NodeType == 6; % end-user water nodes and 23 and 12kv substations
demand_nodeLogic = NodeTable.NodeType == 3; 
demand_nodeIDs = NodeTable.NodeID(demand_nodeLogic); %16-49

num_nodes = size(nodes_can_Fail,1);
num_links = size(LinkTable,1);

Ist_power=[0 1 2 4 6 8];
mean_func_node_Interaction=zeros(length(Ist_power),num_demandNodes);

for i=1:length(Ist_power)
    i=2
    interLinkPowerSupply=148:155; % supply water to power stations
    a_interLinkPowerSupply=a_powerlink.*1;
    RR_interLinkPowerSupply=a_interLinkPowerSupply.*K_powerLink.*(LinkTable.StartNodePGV_ins(interLinkPowerSupply)+...
        LinkTable.EndNodePGV_ins(interLinkPowerSupply))./2;

    LinkTable.LinkPf(interLinkPowerSupply) = 1-exp(-RR_interLinkPowerSupply.*LinkTable.TotLength(interLinkPowerSupply)./1000);
    LinkTable.LinkPf(148:155)=1;  % disable water to power link
    % interdependent links: power nodes demand water supply. 156:164

    interLinkWaterSupply=156:164; % supply power to pumping stations
    a_interLinkWaterSupply=a_waterlink.*Ist_power(i);
    RR_interLinkWaterSupply=a_interLinkWaterSupply.*K_powerLink.*(LinkTable.StartNodePGV_ins(interLinkWaterSupply)+...
        LinkTable.EndNodePGV_ins(interLinkWaterSupply))./2;

    LinkTable.LinkPf(interLinkWaterSupply) = 1-exp(-RR_interLinkWaterSupply.*LinkTable.TotLength(interLinkWaterSupply)./1000);

      
    
    U_nodes = rand(num_simu, num_nodes);
    U_links = rand(num_simu, num_links);% generate num_scenarios*num_links vectors of random #

    probFailure_nodes = repmat(NodeTable.NodePf(nodes_can_Fail_IDs)',num_simu,1);
    probFailure_links = repmat(LinkTable.LinkPf',num_simu,1);
   
    
    % if U < pf_fail, then the node/link is disabled 
    f_nodes = U_nodes <=probFailure_nodes;  % acceptance-rejection sampling to generate pf*num_scenarios
    f_links = U_links <=probFailure_links;
       

    supply_Logic = NodeTable.NodeType == 4 | NodeTable.NodeType == 5 | NodeTable.NodeType == 6; % one directional input dependence. water pumps depend on eletricity
    supply_nodeIDs = NodeTable.NodeID(supply_Logic);
    functionality_xi = zeros(num_simu, size(demand_nodeIDs,1));

    % eletricNodeDemandWaterID=201:208;

          for i_scenario = 1:num_simu   % control+i and shift+tab
             % links first
             links_scen_orig = [s_base t_base];        % original links
             link_fail_scen_i = f_links(i_scenario,:); % identify failed links
             links_scen_i = links_scen_orig(~link_fail_scen_i',:); % remove links that failed in a simulation

             % then the nodes
             nodes_scen_i_orig = nodes_can_Fail';       
             node_fail_scen_i = f_nodes(i_scenario,:);
             nodes_scen_i = nodes_scen_i_orig(node_fail_scen_i);

             % remove any instance of that node from edge list.(Remove failed nodes???)
             remove_firstcolumn_logical = ismember(links_scen_i(:,1) , nodes_scen_i);
             links_scen_1 = links_scen_i(~remove_firstcolumn_logical,:);
             remove_secondcolumn_logical = ismember(links_scen_1(:,2) , nodes_scen_i);
             links_scen_final = links_scen_1(~remove_secondcolumn_logical,:);

             % check connectivity
             G_i = graph();
        %    G_i = addnode(G_i, size(nodeTable.NodeID,1));
             G_i = addnode(G_i, max(NodeTable.NodeID));
             G_i = addedge(G_i, links_scen_final(:,1), links_scen_final(:,2));

             d_i = distances(G_i);

             d_st=d_i(supply_nodeIDs, demand_nodeIDs); %pick out the distances between pairs of storage tanks and distribution nodes.

             num_isNotInf = sum(double(~isinf(d_st)),1);
             % count the number of not inf in each column sum(,1)
             % converts symbolic values to double precision
             functionality_xi(i_scenario,:) = double(~num_isNotInf == 0);   % takes on 1 or 0. 1 is counted as functionable.
             disp(i_scenario)
         end

    mean_func_node_Interaction(i,:) = sum(functionality_xi(:,:),1)/num_simu;
    
end

NodeTable = waterPowerNodeTable;
LinkTable = waterPowerLinkTable;

num_demandNodes=34;

% num_simu= 5000; 
% create base case u-v pairs: the ID of starting and end nodes respectively
s_base = LinkTable.NodeStartID(:);
t_base = LinkTable.NodeEndID(:);

nodes_can_Fail_IDs = NodeTable.NodeType ~= 3;
nodes_can_Fail = NodeTable.NodeID(nodes_can_Fail_IDs);

% demand_nodeLogic = nodeTable.NodeType == 3 | nodeTable.NodeType == 5| nodeTable.NodeType == 6; % end-user water nodes and 23 and 12kv substations
demand_nodeLogic = NodeTable.NodeType == 3; 
demand_nodeIDs = NodeTable.NodeID(demand_nodeLogic); %16-49

num_nodes = size(nodes_can_Fail,1);
num_links = size(LinkTable,1);

mean_func_node_noInteraction=zeros(length(Ist_power),num_demandNodes);

for i=1:length(Ist_power)
    interLinkPowerSupply=148:155;
    a_interLinkPowerSupply=a_powerlink*1;
    RR_interLinkPowerSupply=a_interLinkPowerSupply.*K_powerLink.*(LinkTable.StartNodePGV_ins(interLinkPowerSupply)+...
        LinkTable.EndNodePGV_ins(interLinkPowerSupply))./2;

    LinkTable.LinkPf(interLinkPowerSupply) = 1-exp(-RR_interLinkPowerSupply.*LinkTable.TotLength(interLinkPowerSupply)./1000);

    % interdependent links: power nodes demand water supply. 156:164

    interLinkWaterSupply=156:164;
    a_interLinkWaterSupply=a_waterlink*Ist_power(i);
    RR_interLinkWaterSupply=a_interLinkWaterSupply.*K_powerLink.*(LinkTable.StartNodePGV_ins(interLinkWaterSupply)+...
        LinkTable.EndNodePGV_ins(interLinkWaterSupply))./2;

    LinkTable.LinkPf(interLinkWaterSupply) = 1-exp(-RR_interLinkWaterSupply.*LinkTable.TotLength(interLinkWaterSupply)./1000);

    %-------------------------------------------------
    % no interaction
    LinkTable.LinkPf(148:155)=1;  % disable water to power link
    LinkTable.LinkPf(156:164)=1;  % disable power to water link
    NodeTable.NodePf(50:108)=1;
    
    U_nodes = rand(num_simu, num_nodes);   %
    U_links = rand(num_simu, num_links);   % generate num_scenarios*num_links vectors of random #

    probFailure_nodes = repmat(NodeTable.NodePf(nodes_can_Fail_IDs)',num_simu,1);
    probFailure_links = repmat(LinkTable.LinkPf',num_simu,1);

    % if U < pf_fail, then the node/link is disabled 
    f_nodes = U_nodes <= probFailure_nodes;  % acceptance-rejection sampling to generate pf*num_scenarios
    f_links = U_links <= probFailure_links;
       
    supply_Logic = NodeTable.NodeType == 2;
    supply_nodeIDs = NodeTable.NodeID(supply_Logic);
    
    functionality_xi = zeros(num_simu, size(demand_nodeIDs,1));
    %-----------------------
    % a node in water distribution network is connected to an electrc power
    % node after removing the dependence link: electric depends on water links

    % eletricNodeDemandWaterID=201:208;
         for i_scenario = 1:num_simu   % control+i and shift+tab
             % links first
             links_scen_orig = [s_base t_base];        % original links
             link_fail_scen_i = f_links(i_scenario,:); % identify failed links
             links_scen_i = links_scen_orig(~link_fail_scen_i',:); % remove links that failed in a simulation

             % then the nodes
             nodes_scen_i_orig = nodes_can_Fail';       
             node_fail_scen_i = f_nodes(i_scenario,:);
             nodes_scen_i = nodes_scen_i_orig(node_fail_scen_i);

             % remove any instance of that node from edge list.(Remove failed nodes???)
             remove_firstcolumn_logical = ismember(links_scen_i(:,1) , nodes_scen_i);
             links_scen_1 = links_scen_i(~remove_firstcolumn_logical,:);
             remove_secondcolumn_logical = ismember(links_scen_1(:,2) , nodes_scen_i);
             links_scen_final = links_scen_1(~remove_secondcolumn_logical,:);

             % check connectivity
             G_i = graph();
        %    G_i = addnode(G_i, size(nodeTable.NodeID,1));
             G_i = addnode(G_i, max(NodeTable.NodeID));
             G_i = addedge(G_i, links_scen_final(:,1), links_scen_final(:,2));

             d_i = distances(G_i);

             d_st=d_i(supply_nodeIDs, demand_nodeIDs); %pick out the distances between pairs of storage tanks and distribution nodes.

             num_isNotInf = sum(double(~isinf(d_st)),1);
             % count the number of not inf in each column sum(,1)
             % converts symbolic values to double precision
             functionality_xi(i_scenario,:) = double(~num_isNotInf == 0);   % takes on 1 or 0. 1 is counted as functionable.
             disp(i_scenario)
         end

    mean_func_node_noInteraction(i,:) = sum(functionality_xi(:,:),1)/num_simu;

end

diff=mean_func_node_noInteraction-mean_func_node_Interaction;

% set(0,'defaultaxeslinestyleorder',{'+','o','*','.','x','s','d','^'});

% % mrk=['ko-' 'gs-' 'b*' 'c+' 'm^' 'rv'];
% for i=1:(length(Ist_power)-1)
%     plot(demand_nodeIDs, diff(i,:))
% end
% legend('No Interaction','Ist0','2*Ist0','4*Ist0','6*Ist0','location','Best');

hold on
% plot(demand_nodeIDs, diff(1,:),'ko-')
plot(demand_nodeIDs, diff(2,:),'bs-')
plot(demand_nodeIDs, diff(3,:),'r*--')
plot(demand_nodeIDs, diff(4,:),'kv-')
plot(demand_nodeIDs, diff(5,:),'bp-')
plot(demand_nodeIDs, diff(6,:),'r+-')
legend('Ist0','2*Ist0','4*Ist0','6*Ist0','location','Best');


ax = gca;
ax.YGrid = 'on';
ax.GridLineStyle = '-';

ylim([0 0.6]);
xlabel('Intermediate water delivery node','fontweight','bold','fontsize',12)
ylabel('Difference between mean functionality ratios','fontweight','bold','fontsize',12)
xticks(demand_nodeIDs);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 13)

box off;
ax2 = axes('Position',get(gca,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k');
set(ax2,'YTick', []);
set(ax2,'XTick', []);
box on;
% legend([stem_without,stem_with],{'Without physical interdependence','With physical interdependence'},'location','Best');


%_____________________
% functionality ratio vs link length
hold off

plot(LinkTable.TotLength(demand_nodeIDs),'r+-')


% difference VS dependence strength



% hold on
% stem(demand_nodeIDs, diff(1,:))
% stem(demand_nodeIDs, diff(2,:))


% serviceability ratio

%     % histogram
%     fig_service_histogram=figure('name','Histogram of the Mean Serviceability under Earthquake');
%     hist(mean_service)
%     xlabel('Mean Serviceability')
%     ylabel('Frequncy')

    % pdf
    fig_service_PDF=figure('name','PDF of the Mean Serviceability udner Earthquake');
    [pdf_meanService,pdf_meanService_xi] = ksdensity(mean_service);
    plot(pdf_meanService_xi,pdf_meanService,'LineWidth',3)
    xlabel('Mean serviceability ratio','fontweight','bold','fontsize',12)
    ylabel('Probability density','fontweight','bold','fontsize',12)
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 11)
    % set(gcf, 'Color', 'None') % make the background transparent
    grid on;

    box off;
    ax2 = axes('Position',get(gca,'Position'),...
               'XAxisLocation','top',...
               'YAxisLocation','right',...
               'Color','none',...
               'XColor','k','YColor','k');
    set(ax2,'YTick', []);
    set(ax2,'XTick', []);
    box on;
%     saveas(gcf,fig_service_PDF,'epsc')


% Untitled 5

    i=2
    interLinkPowerSupply=148:155; % supply water to power stations
    a_interLinkPowerSupply=a_powerlink.*1;
    RR_interLinkPowerSupply=a_interLinkPowerSupply.*K_powerLink.*(LinkTable.StartNodePGV_ins(interLinkPowerSupply)+...
        LinkTable.EndNodePGV_ins(interLinkPowerSupply))./2;

    LinkTable.LinkPf(interLinkPowerSupply) = 1-exp(-RR_interLinkPowerSupply.*LinkTable.TotLength(interLinkPowerSupply)./1000);
    LinkTable.LinkPf(148:155)=1;  % disable water to power link
    % interdependent links: power nodes demand water supply. 156:164

    interLinkWaterSupply=156:164; % supply power to pumping stations
    a_interLinkWaterSupply=a_waterlink.*Ist_power(i);
%     RR_interLinkWaterSupply=a_interLinkWaterSupply.*K_powerLink.*(LinkTable.StartNodePGV_ins(interLinkWaterSupply)+...
%         LinkTable.EndNodePGV_ins(interLinkWaterSupply))./2;

RR_interLinkWaterSupply=[0.01 0.02 0.04 0.06];
hold on
Ist=zeros(length(RR_interLinkWaterSupply),length(interLinkWaterSupply));
for i=1:length(RR_interLinkWaterSupply)
    Ist(i,:) = 1-exp(-RR_interLinkWaterSupply(i).*LinkTable.TotLength(interLinkWaterSupply)./1000);
   
    
%     Ist=LinkTable.LinkPf(interLinkWaterSupply);

%     plot(sort(LinkTable.TotLength(interLinkWaterSupply)/1000),sort(LinkTable.LinkPf(interLinkWaterSupply)),'-d')
end
 
% 
hold on
% plot(sort(LinkTable.TotLength(interLinkWaterSupply/1000),sort(LinkTable.LinkPf(interLinkWaterSupply)),'ko-')
plot(sort(LinkTable.TotLength(interLinkWaterSupply)./1000), sort(Ist(1,:)),'ks-');
plot(sort(LinkTable.TotLength(interLinkWaterSupply)./1000), sort(Ist(2,:)),'r*-');
plot(sort(LinkTable.TotLength(interLinkWaterSupply)./1000), sort(Ist(3,:)),'bv-');
plot(sort(LinkTable.TotLength(interLinkWaterSupply)./1000), sort(Ist(4,:)),'mv-');
legend('I0','2*I0','4*I0','6*I0','location','Best');


hold on
% plot(sort(LinkTable.TotLength(interLinkWaterSupply/1000),sort(LinkTable.LinkPf(interLinkWaterSupply)),'ko-')
plot(LinkTable.TotLength(interLinkWaterSupply)./1000, Ist(1,:),'ks-');
plot(LinkTable.TotLength(interLinkWaterSupply)./1000, Ist(2,:),'r*-');
plot(LinkTable.TotLength(interLinkWaterSupply)./1000, Ist(3,:),'bv-');
plot(LinkTable.TotLength(interLinkWaterSupply)./1000, Ist(4,:),'mv-');
legend('I0','2*I0','4*I0','6*I0','location','Best');

ax = gca;
ax.YGrid = 'on';
ax.GridLineStyle = '-';

ylim([0 0.6]);
xlabel('Distance (1000 feet)','fontweight','bold','fontsize',12)
ylabel('Interdependence strength','fontweight','bold','fontsize',12)

box off;
ax2 = axes('Position',get(gca,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k');
set(ax2,'YTick', []);
set(ax2,'XTick', []);
box on;
grid on
% legend([stem_without,stem_with],{'Without physical interdependence','With physical interdependence'},'location','Best');


