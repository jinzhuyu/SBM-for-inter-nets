function ImportanceRemainCompsSort =rank_components_importantce_dynamic(RemainDamagedCompsArray)
    global LinksNotFailedFinalTemp ResilTotWithoutNewInterLinks NodesFailedScenILabelTemp LinksNotFailedScenIMayHaveFailedNodes
    nRemainDamagedComp=length(RemainDamagedCompsArray(:,1));
    ImportanceRemainComps=zeros(nRemainDamagedComp,1);
    ResilRemainComps=zeros(nRemainDamagedComp,1);
    ResilWater=zeros(nRemainDamagedComp,1);
    ResilPower=zeros(nRemainDamagedComp,1);
    
    for iRemainDamagedComp=1:nRemainDamagedComp
        %from previous step. Same for all remain damaged comp this step
        % they are not updated by ranking but after a comp is restored.
        LinksNotFailedFinalTemp1=LinksNotFailedFinalTemp;
        NodesFailedScenILabelTemp1=NodesFailedScenILabelTemp;
        if RemainDamagedCompsArray(iRemainDamagedComp,1)~=RemainDamagedCompsArray(iRemainDamagedComp,2)
            LinksNotFailedFinalTemp1=vertcat(LinksNotFailedFinalTemp1(:,:),RemainDamagedCompsArray(iRemainDamagedComp,1:2));
            [ResilTotTemp,ResilWaterTemp,ResilPowerTemp]=calculate_resi_given_links(LinksNotFailedFinalTemp1);
            ImportanceRemainComps(iRemainDamagedComp)=(ResilTotTemp-ResilTotWithoutNewInterLinks)/ResilTotWithoutNewInterLinks;
            ResilRemainComps(iRemainDamagedComp)=ResilTotTemp;
            ResilWater(iRemainDamagedComp)=ResilWaterTemp;
            ResilPower(iRemainDamagedComp)=ResilPowerTemp;
        else
            % remove the current node from failed nodes
            NodesFailedScenILabelTemp2 = NodesFailedScenILabelTemp1(NodesFailedScenILabelTemp1~=RemainDamagedCompsArray(iRemainDamagedComp,1)); 

            LinksStartNodesFailedLogic = ismember(LinksNotFailedScenIMayHaveFailedNodes(:,1),NodesFailedScenILabelTemp2);
            LinksFailedEndNodesLogic = ismember(LinksNotFailedScenIMayHaveFailedNodes(:,2),NodesFailedScenILabelTemp2);
            LinksWithFailedNodesLogic=LinksStartNodesFailedLogic+LinksFailedEndNodesLogic;
    %       NodesFailedScenILabelTemp=NodesFailedScenILabelTemp1;

            LinksNotFailedFinalTemp2 =  LinksNotFailedScenIMayHaveFailedNodes(~LinksWithFailedNodesLogic,:);

            LinksNotFailedFinalTemp3=union(LinksNotFailedFinalTemp2,LinksNotFailedFinalTemp1,'rows');

            % check connectivity
            [ResilTotTemp,ResilWaterTemp,ResilPowerTemp] = calculate_resi_given_links(LinksNotFailedFinalTemp3);

            ImportanceRemainComps(iRemainDamagedComp)=(ResilTotTemp-ResilTotWithoutNewInterLinks)/ResilTotWithoutNewInterLinks;
            ResilRemainComps(iRemainDamagedComp)=ResilTotTemp;
            ResilWater(iRemainDamagedComp)=ResilWaterTemp;
            ResilPower(iRemainDamagedComp)=ResilPowerTemp;
            %  Resil(iStep)=ResilTotTemp;    
            %  Perhaps in calculating the importance, use Resil without new interlinks.
        end  

    end
    [ImportanceCompsDescend,ImportanceCompsRowIndex]=sort(ImportanceRemainComps,'descend');
    ImportanceRemainCompsSort=[RemainDamagedCompsArray(ImportanceCompsRowIndex,1:2)...
            ResilRemainComps(ImportanceCompsRowIndex) ImportanceCompsDescend...
            ResilWater(ImportanceCompsRowIndex) ResilPower(ImportanceCompsRowIndex)];

end

% rank the components

% select the 1st comp in the ranking

% restore the comp and reconstructure the network
% calculate Resil

% when debugging, walk through one loop step by step.

% keep an eye on the temp variables

% remove the restored componenet, and add it to the importance components
% table

% rank the remaining components

% iterate

