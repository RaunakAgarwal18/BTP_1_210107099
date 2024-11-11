function [X,Moth_fitness] = initialization(SearchAgents_no,dim,lb,ub,fobj,caseNo)

for i=1:dim 
    ub_i=ub(i); 
    lb_i=lb(i);
    X(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;
end

for i2 = 1:SearchAgents_no
    [X(i2,:) , Moth_fitness(i2,:)] = fobj(X(i2,:),caseNo);
end