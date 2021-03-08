function [indexAeqB] = bitwise_equal(A,B)
% return index of elements of vector A that are equal to elements in vector B
%   Detailed explanation goes here
    lenB = length(B);
    indexAeqB = zeros(1,lenB);
    for i=1:lenB
        indexAeqB(i)=find(A==B(i));
    end

end

