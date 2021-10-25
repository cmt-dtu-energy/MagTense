function [TheThreshold] = GetThreshold(TheTensor,FractionOfEntries)
% Compute the threshold corresponding to removing from the tensor
% a given fraction of the non-zero entries

TheFields = fieldnames(TheTensor) ;
for n=1:numel(TheFields)
    ThisField = getfield(TheTensor,TheFields{n}) ;
    if isequal(class(ThisField),'cell')
        Nentries(n) = numel(ThisField{1}) ; 
    else
        Nentries(n) = numel(ThisField) ;
    end
end
NNentries = sum(Nentries) ;
AllEntries = zeros(NNentries,1) ;
for n=1:numel(TheFields)
    ThisField = getfield(TheTensor,TheFields{n}) ;
    if isequal(class(ThisField),'cell')
        AllEntries(sum(Nentries(1:n-1))+1:sum(Nentries(1:n)),1) = reshape(ThisField{1},Nentries(n),1) ;
    else
        AllEntries(sum(Nentries(1:n-1))+1:sum(Nentries(1:n)),1) = reshape(ThisField,Nentries(n),1) ;
    end
end

AllEntries = sort(abs(AllEntries)) ;
AllEntries(AllEntries(:) == 0) =  [] ;
x = [1:numel(AllEntries)]./numel(AllEntries) ;

imax = round(FractionOfEntries*numel(AllEntries)) ;
if imax > numel(AllEntries)
    disp('Threshold fraction >= 1 selected. Turning off demagnetization.')
    TheThreshold = inf;
elseif imax > 0
    TheThreshold = AllEntries(imax) ;
else
    TheThreshold = 0 ;
end
