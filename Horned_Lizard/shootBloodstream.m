%  Programmer: Hernan Peraza    hperaza@ipn.mx
%****************************************************
function o= shootBloodstream(Xbest, X, Max_iter,t)
%------------------------
g= 0.009807;  % 9.807 m/s2   a kilometros    =>  0.009807 km/s2
epsilon= 1E-6;
Vo=  1;%1E-2;
Alpha= pi/2;
%-------------------------
 o= ( Vo * cos(Alpha*t/Max_iter)+epsilon) * Xbest + (Vo * sin(Alpha-Alpha*t/Max_iter)-g+epsilon) * X;
end
%**********************************************************************

