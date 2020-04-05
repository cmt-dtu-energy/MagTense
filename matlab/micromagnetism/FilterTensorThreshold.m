function [TheTensor] = FilterTensorThreshold(TheTensor,TheThreshold)
% remove (i.e. set to zero) from the tensor the entries whose absolute
% value is smaller than the given threshold.

TheFields = fieldnames(TheTensor) ;

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