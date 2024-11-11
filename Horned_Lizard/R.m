%  Programmer: Hernan Peraza    hperaza@ipn.mx
%****************************************************
function [r1,r2,r3,r4]= R(NP)
     band=1;
        while band
            r1= round(1+ (NP-1)*rand());
            r2= round(1+ (NP-1)*rand());
            r3= round(1+ (NP-1)*rand()); 
            r4= round(1+ (NP-1)*rand()); 
            if (r1 ~= r2) && ( r2 ~= r3)  && ( r1 ~= r3) && ( r4 ~= r3) && ( r4 ~= r2)&& ( r1 ~= r4)
                band= 0;
            end
        end    
end
%******************************************************

