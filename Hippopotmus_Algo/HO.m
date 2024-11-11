%% Designed and Developed by Mohammad Hussien Amiri and Nastaran Mehrabi Hashjin
function[Best_Fitness_HO,Best_Hippo_HO,Convergence_curve_HO,NumberOfEval] = HO(Hippo_Pos,Hippo_fitness,SearchAgents,Max_iterations,lowerbound,upperbound,dimension,fitness,RandomAgentsRequired,caseNo,maxEval)
 
    if RandomAgentsRequired > 0
        [Hippo_Pos,Hippo_fitness] = initialization(SearchAgents,dimension,lowerbound,upperbound,fitness,caseNo);   % Initializing Search Agents if RandomAgentsRequired >=1
    end
    NumberOfEval = SearchAgents;
    Convergence_curve_HO = inf*ones(1,Max_iterations+1);
    Best_Hippo_HO        = zeros(1,dimension);
    Best_Fitness_HO      = 0;
    Convergence_curve_HO(1) = min(Hippo_fitness);
    %% Main Loop
    
    for t=1:Max_iterations

        %% Update the Best Condidate Solution/ Dominant hippo
        [best,location]=min(Hippo_fitness);
        if t==1
            Xbest = Hippo_Pos(location,:);                                  % Optimal location
            fbest = best;                                                   % The optimization objective function
        elseif best<fbest
            fbest = best;
            Xbest = Hippo_Pos(location,:);
        end
        
          % Half population goes through phase 1 and other half through
          % phase 2
        for i=1:SearchAgents/2
    
            %% Phase1: The hippopotamuses position update in the river or pond (Exploration) (immature mal and female hippo)
            Dominant_hippopotamus = Xbest;
            I1                    = randi([1,2],1,1);                       % Integer : 1 or 2
            I2                    = randi([1,2],1,1);                       % Integer : 1 or 2
            Ip1                   = randi([0,1],1,2);                       % Ip1(1) and Ip1(2) either 0 or 1
            RandGroupNumber       = randperm(SearchAgents,1);               % selects a random number between 1 and number of search agents
            RandGroup             = randperm(SearchAgents,RandGroupNumber); % vector containing RandGroupNumber unique random indices.
    
            % Mean of Random Group
            MeanGroup = mean(Hippo_Pos(RandGroup,:)).*(length(RandGroup)~=1) + Hippo_Pos(RandGroup(1,1),:)*(isscalar(RandGroup)); % <- Randgroup has more than 1 hippo
                             % case where RandGroup has only 1 hippo
            y1 = rand(1,1);
            Alfa{1,:} = (I2*rand(1,dimension)+(~Ip1(1)));
            Alfa{2,:} = 2*rand(1,dimension)-1;
            Alfa{3,:} = rand(1,dimension);
            Alfa{4,:} = (I1*rand(1,dimension)+(~Ip1(2)));
            Alfa{5,:} = rand;
            A         = Alfa{randi([1,5],1,1),:};  % random value from the alfa vector
            B         = Alfa{randi([1,5],1,1),:};  % random value from the alfa vector
            X_P1(i,:) = Hippo_Pos(i,:)+y1.*(Dominant_hippopotamus-I1.*Hippo_Pos(i,:));  % X_P1 represents the male hippo.
            T         = exp(-t/Max_iterations);
            if T>0.6  % immature hippo moved away from his mother
                X_P2(i,:) = Hippo_Pos(i,:)+A.*(Dominant_hippopotamus-I2.*MeanGroup);           % x_P2 represents the female/immature hippo.
            else
                if rand()>0.5  % immature hippo is away from mother but still near the herd
                    X_P2(i,:) = Hippo_Pos(i,:)+B.*(MeanGroup-Dominant_hippopotamus);
                else           % immature hippo moved away from the herd
                    X_P2(i,:) = ((upperbound-lowerbound)*rand+lowerbound);
                end
            end

            X_P2(i,:) = min(max(X_P2(i,:),lowerbound),upperbound);
            L         = X_P1(i,:);
            [L,F_P1(i)]   = fitness(L,caseNo);
            NumberOfEval = NumberOfEval+1;
            if(F_P1(i)<Hippo_fitness(i))
                Hippo_Pos(i,:)   = X_P1(i,:);
                Hippo_fitness(i) = F_P1(i);
            end
    
            L2      = X_P2(i,:);
            [L2,F_P2(i)] = fitness(L2,caseNo);
            NumberOfEval = NumberOfEval+1;
            if(F_P2(i)<Hippo_fitness(i))
                Hippo_Pos(i,:)   = X_P2(i,:);
                Hippo_fitness(i) = F_P2(i);
            end
        end

        %% Phase 2: Hippopotamus defense against predators (Exploration) 
        for i = 1+SearchAgents/2:SearchAgents
            predator        = lowerbound+rand(1,dimension).*(upperbound-lowerbound);  %randomly generated predator
            L               = predator;
            [L,F_HL]        = fitness(L,caseNo);
            NumberOfEval    = NumberOfEval+1;
            distance2Leader = abs(predator-Hippo_Pos(i,:));  %distance of predator from the leader of herd
            
            % each line below generates a uniform random number
            b               = unifrnd(2,4,[1 1]);     % random number between 2 and 4 with equal probability
            c               = unifrnd(1,1.5,[1 1]);   % random number between 1 and 1.5 with equal probability
            d               = unifrnd(2,3,[1 1]);     % random number between 2 and 3 with equal probability
            l               = unifrnd(-2*pi,2*pi,[1 1]); % random number between -2pi and 2pi with equal probability
            
            RL              = 0.05*levy(SearchAgents,dimension,1.5);  % random vector based on levy distribution representing levi movement.
    
            if Hippo_fitness(i)>F_HL   % Predator is far from hippo
                X_P3(i,:) = RL(i,:).*predator+(b./(c-d*cos(l))).*(1./distance2Leader);
            else                       % Predator is near the hippo
                X_P3(i,:) = RL(i,:).*predator+(b./(c-d*cos(l))).*(1./(2.*distance2Leader+rand(1,dimension)));
            end

            X_P3(i,:) = min(max(X_P3(i,:),lowerbound),upperbound);
            L         = X_P3(i,:);
            [L,F_P3(i)]   = fitness(L,caseNo);
            
            NumberOfEval = NumberOfEval+1;
            if(F_P3(i)<Hippo_fitness(i))
                Hippo_Pos(i,:) = X_P3(i,:);
                Hippo_fitness(i) = F_P3(i);
            end
        end
    
    
        %% Phase 3: Hippopotamus Escaping from the Predator (Exploitation)
        for i = 1:SearchAgents
            LO_LOCAL  = (lowerbound./t);
            HI_LOCAL  = (upperbound./t);
            Alfa{1,:} = 2*rand(1,dimension)-1;
            Alfa{2,:} = rand(1,1);
            Alfa{3,:} = randn;  % random value from a standard normal distribution (mean 0, standard deviation 1)
            D         = Alfa{randi([1,3],1,1),:}; % selects one of the three cells in Alfa to assign to D.
            X_P4(i,:) = Hippo_Pos(i,:)+(rand(1,1)).*(LO_LOCAL+D.* (HI_LOCAL-LO_LOCAL));  % Generating a nearby safe position
            X_P4(i,:) = min(max(X_P4(i,:),lowerbound),upperbound);
            L         = X_P4(i,:);
            [L,F_P4(i)]   = fitness(L,caseNo);
            NumberOfEval = NumberOfEval+1;
            if(F_P4(i)<Hippo_fitness(i))
                Hippo_Pos(i,:) = X_P4(i,:);
                Hippo_fitness(i) = F_P4(i);
            end
        end
           
        % Find the global best solution
        [Best_Score, location]    = min(Hippo_fitness);    
        Best_Hippo_HO             = Hippo_Pos(location,:);
        Best_Fitness_HO           = Best_Score;
        Convergence_curve_HO(t+1) = Best_Score;
    end
end

