%  Programmer: Hernan Peraza    hperaza@ipn.mx
%****************************************************
function o=checkBoundaries(X, lb, ub)
  ubx=[];
  lbx=[];
  n= size(X,2);    
  bound= size(lb,2);
   if bound==1
      for i=1:n
         ubx(i)= ub;   
         lbx(i)= lb;
      end    
   end
   if bound >1
      for i=1:n
        ubx(i)= ub(i);   
        lbx(i)= lb(i);
      end
   end
 X(X>ub)= ubx(X>ub);
 X(X<lb)= lbx(X<lb);
 o=X;
end
%********************************************
