%  Programmer: Hernan Peraza    hperaza@ipn.mx
%****************************************************

% color        Hexadecimal  Decimal   normalizado
% Darkening1    F9D100     16371968      0.1
% Darkening2    9398BF      9672895   0.43767961
% Lightening1   7636608     747080    0.54090919
% Lightening2   161617      1447447    0.9
%***************************************************
% color	Hexadecimal	Decimal	normalizado
% Lightening 1	e8e8e8	15263976	0.0000000
% Lightening 2	9398BF	9672895	0.4046661
% Darkening 1	763660	7747080	0.5440510
% Darkening 2	161617	1447447	1.0000000
%*************************************************************************
function o =  Skin_darkening_or_lightening(Xbest,  X, SearchAgents_no)
 darkening= [0.0, 0.4046661];
 lightening= [0.5440510 , 1.0];

    dark1= darkening(1) + (darkening (2) - darkening(1))* rand();
    dark2= darkening(1) + (darkening (2) - darkening(1))* rand();
    light1= lightening(1) + (lightening(2)-lightening(1))* rand();
    light2= lightening(1) + (lightening(2)-lightening(1))* rand();

  [r1, r2, r3, r4]= R(SearchAgents_no);

    if(getBinary)
         o= Xbest +   light1*sin((X(r1,:)-X(r2,:))/2) - ((-1)^getBinary) * light2*sin((X(r3,:)-X(r4,:))/2);
    else
         o= Xbest +    dark1*sin((X(r1,:)-X(r2,:))/2) - ((-1)^getBinary) * dark2*sin((X(r3,:)-X(r4,:))/2);
    end
end
%*************************************************************************

