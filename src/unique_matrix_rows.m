function A=unique_matrix_rows(A)


[UniA,Index]=unique(A,'rows');

DupIndex=setdiff(1:size(A,1),Index);
A(DupIndex,:)=[];
return 