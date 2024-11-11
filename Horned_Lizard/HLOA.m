function [vMin,theBestVct,Convergence_curve,NumberOfEval]=HLOA(Positions,Fitness,SearchAgents_no,Max_iter,lb,ub,dim,fobj,RandonRequired,caseNo,maxEval)

% Positions=initialization_HLOA(SearchAgents_no,dim,ub,lb); % Position - searchAgentNo. * dim (no. of decision variables)
if(RandonRequired > 0)
    Fitness = zeros(SearchAgents_no,1);
    for i=1:size(Positions,1)
      [Positions(i,:),Fitness(i)] = fobj(Positions(i,:),caseNo); % get fitness  
    end
end
NumberOfEval = SearchAgents_no;
[vMin,minIdx] = min(Fitness);  % the min fitness value vMin and the position minIdx
theBestVct = Positions(minIdx,:);  % the best vector 
[vMax,maxIdx] = max(Fitness); % the max fitness value vMax and the position maxIdx
Convergence_curve = inf*ones(1,Max_iter+1);
Convergence_curve(1) = vMin;
alphaMelanophore = alpha_melanophore(Fitness,vMin,vMax);
% Main 
v = zeros(SearchAgents_no,dim);
for t=1:Max_iter  
   for r=1:SearchAgents_no
     if( 0.5 < rand)   % se mimetiza 
             v(r,:)= mimicry(theBestVct, Positions, Max_iter, SearchAgents_no, t);
     else
            if(mod(t,2))
                v(r,:)= shootBloodstream(theBestVct, Positions(r,:), Max_iter,t); 
            else
                v(r,:)= randomWalk(theBestVct,Positions(r,:));
            end
     end   
     Positions(maxIdx,:)= Skin_darkening_or_lightening(theBestVct, Positions, SearchAgents_no);
     if (alphaMelanophore(r) <= 0.3)
         v(r,:)= remplaceSearchAgent(theBestVct, Positions,SearchAgents_no);
    end        
     %----------------------------------------------------------------------        
     % Return back the search agents that go beyond the boundaries of the search space
     v(r,:)=checkBoundaries(v(r,:), lb, ub);
     % Evaluate new solutions
     [v(r,:),Fnew] = fobj(v(r,:),caseNo);
     NumberOfEval = NumberOfEval+1;
     % Update if the solution improves
     if Fnew <= Fitness(r)
        Positions(r,:)= v(r,:);
        Fitness(r)= Fnew;
     end
     if Fnew <= vMin
         theBestVct= v(r,:);
         vMin= Fnew;
     end 
   end
   %update max and alpha-melanophore
   [vMax,maxIdx]= max(Fitness);
   alphaMelanophore = alpha_melanophore(Fitness,vMin,vMax);
   Convergence_curve(t+1)= vMin; 
 end
%***********************************[End HLOA Algorithm]
