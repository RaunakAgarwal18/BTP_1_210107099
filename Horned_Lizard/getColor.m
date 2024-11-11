%  Programmer: Hernan Peraza    hperaza@ipn.mx
%****************************************************
function [c1,c2]= getColor(colorPalette)
 band=1;
         while band
            c1 = colorPalette(1,randi([1 30],1,1));
            c2 = colorPalette(1,randi([1 30],1,1)) ;
            if (c1 ~= c2)   
                band= 0;
            end
        end    
end
%*************************************************
