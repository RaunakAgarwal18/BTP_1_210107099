%  Programmer: Hernan Peraza    hperaza@ipn.mx
%****************************************************
function o = remplaceSearchAgent(Xbest, X, SearchAgents_no)
band=1;
  while band
       r1= round(1+ (SearchAgents_no-1)*rand());
       r2= round(1+ (SearchAgents_no-1)*rand());
       if r1 ~= r2
            band= 0;
       end
  end
  o=  Xbest + (X(r1,:)-((-1)^getBinary)*X(r2,:))/2; 
end
%*********************************************************************
