%  Moth Swarm Algorithm (MSA)                                                            
function [Best_Fitness_MSA,Best_Moth_MSA,Convergence_curve_MSA,NumberOfEval]=MSA(Moth_pos,Moth_fitness,SearchAgents,Max_iterations,lowerbound,upperbound,dimension,fitness,RandomAgentsRequired,caseNo)
    
    Nc = SearchAgents/5; 
    % optimal solution is the position of light source
    if RandomAgentsRequired > 0
        [Moth_pos,Moth_fitness] = initialization(SearchAgents,dimension,lowerbound,upperbound,fitness,caseNo);   % Initializing Search Agents if RandomAgentsRequired >=1
    end
    NumberOfEval = SearchAgents;
    [Moth_fitness,location]  = sort(Moth_fitness);                      % Sorting the fitness value in ascending order
    Moth_pos                 = Moth_pos(location,:);                    % Organizing the search Agaents based on sorted fitness values
    Convergence_curve_MSA    = inf*ones(1,Max_iterations+1);
    Best_Moth_MSA            = zeros(1,dimension);
    Best_Fitness_MSA         = 0;
    Convergence_curve_MSA(1) = min(Moth_fitness);
    
    % Main loop
    for t = 1:Max_iterations 
        pMoth_pos     = Moth_pos;
        pMoth_fitness = Moth_fitness; 
        % 3 groups of moths - Pathfinder, prospectors and onlookers
        %____________________ 1. Reconnaissance phase : Pathfinder moths - Search the optimization search space _______________    
        
        Lights        = Moth_pos(1:Nc,:);   
        Light_fitness = Moth_fitness(1:Nc);      % best moths considered to be pathfinders, Moth_fitness is sorted in ascending order
        
        [Lights,Light_fitness] = Pathfinders(Lights,Light_fitness,Nc,fitness,upperbound,lowerbound,caseNo); % improve pathfinders
        NumberOfEval = NumberOfEval + Nc;
        Moth_fitness(1:Nc) = Light_fitness;      
        Moth_pos(1:Nc,:)   = Lights;             % insertion in the swarm
        
        % Sharing of luminescence intensities using Eqs.(28) and (29). 
        for i = 1:Nc
            if Light_fitness(i) >= 0
                Light_fitness(i) = 1/(1+Light_fitness(i));
            else
                Light_fitness(i) = 1+abs(Light_fitness(i));
            end
        end
        
        for i=1:Nc
            probability(i) = Light_fitness(i)/sum(Light_fitness);
        end 
        
        [~,R]     = histc(rand(1,SearchAgents),cumsum([0;probability(:)]));
        a1        = 1:Nc;
        new_Light = Lights(a1(R),:);
        
        %__________________2. Transverse orientation phase : prospector moths_________________________________    
        % inspired from Moth-flame Optimization (MFO) [14]: with each variable as an integrated unit
        
        Predictor_no = round((SearchAgents-Nc)*(1-t/(Max_iterations))); %No. of prospectors using Eq.(30).
        a            = -1-t/Max_iterations;
        
        for i=Nc+1:Predictor_no+Nc
               tt            = (a-1)*rand()+1; 
               spiral        = exp(tt).*cos(tt.*2*pi);    
               D_to_Light    = abs(new_Light(i,:)-Moth_pos(i,:));
               Moth_pos(i,:) = D_to_Light.*spiral+new_Light(i,:);       % updating prospectors using Eq.(31)          
        end
        
        % Return back prospectors that go beyond the boundaries
        Moth_pos(Nc+1:Predictor_no+Nc,:) = Bound_Checking(Moth_pos(Nc+1:Predictor_no+Nc,:),upperbound,lowerbound);
        
        %Fitness evaluation for prospectors
        for i2 = Nc+1:Predictor_no+Nc
            [Moth_pos(i2,:),Moth_fitness(i2,:)] = fitness(Moth_pos(i2,:),caseNo);
        end
        NumberOfEval = NumberOfEval + Predictor_no;
        [~, location] = min(Moth_fitness);    
        gx            = Moth_pos(location,:);
        
        %_______________________ 3. Celestial navigation phase: onlooker moths____________________
        % 
        
        % 3.1 Gaussian walks

         x3 = round((SearchAgents-Predictor_no-Nc)*1/2);
         j2 = SearchAgents;
        
         for j1=1:x3
            j2 = j1+Predictor_no+Nc;
    
            % updating prospectors using Eqs.(32-33)
            Moth_pos(j2,:)   = normrnd(gx(1,:),(log(t)/t)*(abs((Moth_pos(j2,:)-gx(1,:)))), [1 size(Moth_pos,2)])+(randn*gx(1,:)-randn*Moth_pos(j2,:));    %gaussian(large step away)
            out              = Moth_pos(j2,:)<lowerbound | Moth_pos(j2,:) > upperbound;
            Moth_pos(j2,out) = normrnd(pMoth_pos(j2,out),(abs((rand*pMoth_pos(j2,out)-rand*gx(1,out)))), [1 size(Moth_pos(j2,out),2)]);                 %(step around);    
         
         end
           
        % 3.2 Associative learning mechanism with immediate memory (ALIM)
        
        x3 = SearchAgents-j2; % no. of onlookers moves with ALIM
        
        for j1 = 1:x3 
            j3 = j1+j2;
            % updating prospectors using Eq.(34)
            Moth_pos(j3,:) = Moth_pos(j3,:)+0.001*unifrnd(lowerbound-Moth_pos(j3,:),upperbound-Moth_pos(j3,:))+(2*t/Max_iterations)*rand(size(gx)).*(gx(1,:)-Moth_pos(j3,:))+(1-t/Max_iterations)*rand(size(gx)).*(new_Light(j3,:)-Moth_pos(j3,:));%step away like PSO
        end

        % Apply Position Limits for onlookers
        Moth_pos(Predictor_no+Nc+1:SearchAgents,:) = Bound_Checking(Moth_pos(Predictor_no+Nc+1:SearchAgents,:),upperbound,lowerbound);
        
        %Fitness evaluation for onlookers 
        for i2 = Predictor_no+Nc+1:SearchAgents
            [Moth_pos(i2,:),Moth_fitness(i2,:)] = fitness(Moth_pos(i2,:),caseNo);
        end
        NumberOfEval = NumberOfEval + SearchAgents - (Predictor_no + Nc);
        %===================================================================== 
        % select the best moths
        [Moth_fitness,I] = unique([Moth_fitness, pMoth_fitness],'first');  
        Moth_fitness     = Moth_fitness(1:SearchAgents);
        dMoth_pos        = [Moth_pos; pMoth_pos];
        Moth_pos         = dMoth_pos(I(1:SearchAgents),:);
    
        % Find the global best solution
        [Best_score, location]     = min(Moth_fitness);    
        Best_Moth_MSA              = Moth_pos(location,:);
        Best_Fitness_MSA           = Best_score;
        Convergence_curve_MSA(t+1) = Best_score;

    end
end