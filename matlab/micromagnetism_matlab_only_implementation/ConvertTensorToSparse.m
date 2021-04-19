function TheTensor = ConvertTensorToSparse(TheTensor)

TheFields = fieldnames(TheTensor) ;


for n=1:numel(TheFields)
    ThisField = getfield(TheTensor,TheFields{n}) ;
        if isequal(class(ThisField),'cell')
            ThisField{1} = sparse(ThisField{1}) ;
            TheTensor = setfield(TheTensor,TheFields{n},ThisField) ;
        else
            TheTensor = setfield(TheTensor,TheFields{n},sparse(ThisField)) ;
        end
end