function [theString] = DoTheUnitThing(theNumber,theUnit)

if theNumber==0
    theString = ['0 ',theUnit] ;
else
    TheSign = sign(theNumber) ;
    SignSym = {'-','','','',''} ;
    theNumber = abs(theNumber) ;
    thePrefix = {'a','f','p','n','u','m','','k','M','G','T','P','E'} ;
    theLog10 = round(log10(theNumber)/3)*3 ;
    theNumber = theNumber/(10^theLog10) ;
    try
        theString = [SignSym{TheSign+2},num2str(theNumber,'%0.3f'),' ',thePrefix{theLog10/3+7},theUnit] ;
    catch
        theString = [SignSym{TheSign+2},num2str(theNumber,'%0.3f'),'*10^{',num2str(theLog10),'} ',theUnit] ;
    end
end

end