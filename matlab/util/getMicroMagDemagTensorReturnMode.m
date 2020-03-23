function returnMode = getMicroMagDemagTensorReturnMode( mode )

    
    switch( mode )
        case 'donot'%do not return N
            returnMode = int32(1);            
        case 'memory' %attempt to return in memory
            returnMode = int32(2);            
        otherwise %return as a file and the length of the file is equal to mode (converted to an int32)
            returnMode = int32(length(mode));            
    end

end