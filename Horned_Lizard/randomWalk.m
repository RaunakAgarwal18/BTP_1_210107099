%  Programmer: Hernan Peraza    hperaza@ipn.mx
%****************************************************
function o = randomWalk(Xbest,X)
  e=  CauchyRand(0,1);
  walk= -1 + 2 * rand(); % -1 < d < 1
  o= Xbest + walk*(0.5-e)*X ;
end
%*************************************************
