function [Best_Fitness_TLBO,Best_Agent_TLBO,Convergence_curve_TLBO,NumberOfEval] = TLBO(Agent_pos,Agent_fitness,SearchAgents,Max_iterations,lowerbound,upperbound,dimension,fitness,RandomAgentsRequired,caseNo)
    
    if RandomAgentsRequired > 0
        [Agent_pos,Agent_fitness] = initialization(SearchAgents,dimension,lowerbound,upperbound,fitness,caseNo);   % Initializing Search Agents if RandomAgentsRequired >=1
    end
    NumberOfEval = SearchAgents;
    Convergence_curve_TLBO = inf*ones(1,Max_iterations+1);
    Best_Agent_TLBO = zeros(1,dimension);
    Best_Fitness_TLBO = 0;
    %% Starting of TLBO
    Convergence_curve_TLBO(1) = min(Agent_fitness);
    %% Iteration loop
    for t = 1: Max_iterations

        for i = 1:SearchAgents

            %% Teacher Phase
            Xmean   = mean(Agent_pos);                                         % Determining mean of the population
            [~,ind] = min(Agent_fitness);                                      % Detemining the location of the teacher
            Xbest   = Agent_pos(ind,:);                                        % Copying the solution acting as teacher
            TF      = randi([1 2],1,1);                                        % Generating either 1 or 2 randomly for teaching factor
            Xnew    = Agent_pos(i,:) + rand(1,dimension).*(Xbest - TF*Xmean);  % Generating the new solution
            Xnew    = min(upperbound, Xnew);                                   % Bounding the violating variables to their upper bound
            Xnew    = max(lowerbound, Xnew);                                   % Bounding the violating variables to their lower bound
            [Xnew,fnew]    = fitness(Xnew,caseNo);                                           % Evaluating the fitness of the newly generated solution
            NumberOfEval = NumberOfEval+1;
            if (fnew < Agent_fitness(i))            % Greedy selection
                Agent_pos(i,:)   = Xnew;            % Include the new solution in population
                Agent_fitness(i) = fnew;            % Include the fitness function value of the new solution in population
            end

            %% Learner Phase
            p = randi([1 SearchAgents],1,1);        % Selection of random parter

            %% Ensuring that the current member is not the partner
            while i == p
                p = randi([1 SearchAgents],1,1);    % Selection of random parter
            end

            if Agent_fitness(i)< Agent_fitness(p)   % Select the appropriate equation to be used in Learner phase
                Xnew = Agent_pos(i,:) + rand(1, dimension).*(Agent_pos(i,:) - Agent_pos(p,:));  % Generating the new solution
            else
                Xnew = Agent_pos(i,:) - rand(1, dimension).*(Agent_pos(i,:) - Agent_pos(p,:));  % Generating the new solution
            end

            Xnew = min(upperbound, Xnew);       % Bounding the violating variables to their upper bound
            Xnew = max(lowerbound, Xnew);       % Bounding the violating variables to their lower bound
            [Xnew,fnew] = fitness(Xnew,caseNo);        % Evaluating the fitness of the newly generated solution
            NumberOfEval = NumberOfEval+1;
            if (fnew < Agent_fitness(i))        % Greedy selection
                Agent_pos(i,:) = Xnew;          % Include the new solution in population
                Agent_fitness(i) = fnew;        % Include the fitness function value of the new solution in population
            end
        end
        Convergence_curve_TLBO(t+1) = min(Agent_fitness);      % Storing the best value of each iteration
    end
    [Best_Fitness_TLBO,ind] = min(Agent_fitness);
    Best_Agent_TLBO = Agent_pos(ind,:);
end
