% find the population of the cencus tract where each node is located in.

PowerNodeAttr = xlsread('PowerNodesAttributes.xlsx');

CencusPopul = xlsread('cencusPopulation.xlsx');
PowerNodeAttr(:,11)=zeros(length(PowerNodeAttr(:,7)),1);

for ii=1:length(PowerNodeAttr(:,7))
    row_same=find(CencusPopul(:,1)==round(PowerNodeAttr(ii,7),0));
    PowerNodeAttr(ii,11)=CencusPopul(row_same,2);
end


WaterNodeAttr = xlsread('WaterNodeAttributes.xlsx');

WaterNodeAttr(:,8)=zeros(length(WaterNodeAttr(:,5)),1);

for ii=1:length(WaterNodeAttr(:,5))
    row_same=find(CencusPopul(:,1)==round(WaterNodeAttr(ii,5),0));
    WaterNodeAttr(ii,8)=CencusPopul(row_same,2);
end