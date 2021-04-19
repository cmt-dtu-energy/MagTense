function TheTensor = ApplyFourierTransformToTensor(TheTensor,UU,UU1)

TheFields = fieldnames(TheTensor) ;


for n=1:numel(TheFields)
    ThisField = getfield(TheTensor,TheFields{n}) ;
    if isequal(class(ThisField),'cell')
        ThisField{1} = UU*ThisField{1}*UU1 ;
    else
        ThisField = UU*ThisField*UU1 ;
    end
    TheTensor = setfield(TheTensor,TheFields{n},ThisField) ;
    
end