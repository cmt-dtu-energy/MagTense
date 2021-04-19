function [TheTensor,TheRelThreshold] = FilterTensor(TheTensor,FractionOfEntries)

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
x = [1:NNentries]./NNentries ;
% imax = max(find(x<FractionOfEntries)) ;
imax = round(FractionOfEntries*NNentries) ;
if imax > 0
    TheThreshold = AllEntries(imax) ;
    %     figure ;
    % plot([1:NNentries]./NNentries,AllEntries./AllEntries(end),'.-')
    % hold on
    % plot(x(imax),TheThreshold/AllEntries(end),'or')
    
else
    TheThreshold = 0 ;
end
TheRelThreshold = TheThreshold/AllEntries(end) ;
for n=1:numel(TheFields)
    ThisField = getfield(TheTensor,TheFields{n}) ;
    if isequal(class(ThisField),'cell')
        Iout = abs(ThisField{1})<=TheThreshold ;
        ThisField{1}(Iout) = 0 ;
        TheTensor = setfield(TheTensor,TheFields{n},ThisField) ;
    else
        Iout = abs(ThisField)<=TheThreshold ;
        ThisField(Iout) = 0 ;
        TheTensor = setfield(TheTensor,TheFields{n},ThisField) ;
    end
    
end