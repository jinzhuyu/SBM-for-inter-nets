function [dataNorm] = truncated_normalize(data,normLowBound)
% normalized a vector or each column of a datarix to [normLowerBound,1]
% input: a datarix or a vector
    if ismatrix(data)==1
        dataNorm = zeros(size(data));
        for iRow = 1:length(data(:,1))
%             colMax= max(data(iRow,:));
            colMin= min(data(iRow,:));
            colRange = range(data(iRow,:));
            dataNorm(iRow,:) = (data(iRow,:)-colMin+normLowBound*colRange)./...
                              ((1+normLowBound)*colRange);
        end
    elseif isvector(data)==1
%         dataMax= max(data);
        dataMin= min(data);
        dataRange = range(data);
        dataNorm = (data-dataMin+normLowBound*dataRange)./...
                  ((1+normLowBound)*dataRange);        
    else
        error('The input of the truncated normalization function should be a datarix or a vector');
    end
end

