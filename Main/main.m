clc
clear                  
SearchAgents = 50;                        
Max_iterations = 500;
nRuns = 200;
AgentsReq = 1;
cases = 4;

[product,l,m,h,il,im,ih,cl,cm,ch,SP,rm1,rm2,rm3,nProcess] = ProductionPlanningData;
fitness    = @SKS_ProductionPlanningC;
dimension  = length(l);
lowerbound = zeros(1,dimension);
upperbound = h';

%HLOA
Best_Fitness_HLOA_nRuns = NaN(cases,nRuns,1);
Best_Sol_HLOA_nRuns = NaN(cases,nRuns,dimension);
Convergence_curve_HLOA_nRuns = NaN(cases,nRuns,Max_iterations+1);

%HO
Best_Fitness_HO_nRuns = NaN(cases,nRuns,1);
Best_Sol_HO_nRuns = NaN(cases,nRuns,dimension);
Convergence_curve_HO_nRuns = NaN(cases,nRuns,Max_iterations+1);

%MSA
Best_Fitness_MSA_nRuns = NaN(cases,nRuns,1);
Best_Sol_MSA_nRuns = NaN(cases,nRuns,dimension);
Convergence_curve_MSA_nRuns = NaN(cases,nRuns,Max_iterations+1);

%TLBO
Best_Fitness_TLBO_nRuns = NaN(cases,nRuns,1);
Best_Sol_TLBO_nRuns = NaN(cases,nRuns,dimension);
Convergence_curve_TLBO_nRuns = NaN(cases,nRuns,Max_iterations+1);

%MSA_HO_Hybrid
Best_Fitness_MSA_HO_nRuns = NaN(cases,nRuns,1);
Best_Sol_MSA_HO_nRuns = NaN(cases,nRuns,dimension);
Convergence_curve_MSA_HO_nRuns = NaN(cases,nRuns,Max_iterations+1);


for i = 1:cases  
    for runs = 1:nRuns
        [X,Fitness_Values]=initialization(SearchAgents,dimension,lowerbound,upperbound,fitness,i);
        fprintf("Running Case %d Run %d \n",i,runs);
        %MSA_HO_Hybrid
        [Best_Fitness_MSA_HO_Hybrid,Best_Sol_MSA_HO_Hybrid,Convergence_curve_MSA_HO_Hybrid,NumberOfEval_MSA_HO] = MSA_HO_Hybrid(X,Fitness_Values,SearchAgents,Max_iterations,lowerbound,upperbound,dimension,fitness,AgentsReq,i);
        Best_Fitness_MSA_HO_nRuns(i,runs) = Best_Fitness_MSA_HO_Hybrid;
        Best_Sol_MSA_HO_nRuns(i,runs,:) = Best_Sol_MSA_HO_Hybrid;
        Convergence_curve_MSA_HO_nRuns(i,runs,:) = Convergence_curve_MSA_HO_Hybrid;

        %MSA
        [Best_Fitness_MSA,Best_Sol_MSA,Convergence_curve_MSA,NumberOfEval_MSA] = MSA(X,Fitness_Values,SearchAgents,Max_iterations,lowerbound,upperbound,dimension,fitness,AgentsReq,i);
        Best_Fitness_MSA_nRuns(i,runs) = Best_Fitness_MSA;
        Best_Sol_MSA_nRuns(i,runs,:) = Best_Sol_MSA;
        Convergence_curve_MSA_nRuns(i,runs,:) = Convergence_curve_MSA;

        %HO
        [Best_Fitness_HO,Best_Sol_HO,Convergence_curve_HO,NumberOfEval_HO] = HO(X,Fitness_Values,SearchAgents,Max_iterations,lowerbound,upperbound,dimension,fitness,AgentsReq,i);
        Best_Fitness_HO_nRuns(i,runs) = Best_Fitness_HO;
        Best_Sol_HO_nRuns(i,runs,:) = Best_Sol_HO;
        Convergence_curve_HO_nRuns(i,runs,:) = Convergence_curve_HO;

        %TLBO
        [Best_Fitness_TLBO,Best_Sol_TLBO,Convergence_curve_TLBO,NumberOfEval_TLBO] = TLBO(X,Fitness_Values,SearchAgents,Max_iterations,lowerbound,upperbound,dimension,fitness,AgentsReq,i);
        Best_Fitness_TLBO_nRuns(i,runs) = Best_Fitness_TLBO;
        Best_Sol_TLBO_nRuns(i,runs,:) = Best_Sol_TLBO;
        Convergence_curve_TLBO_nRuns(i,runs,:) = Convergence_curve_TLBO;

        %HLOA
        [Best_Fitness_HLOA,Best_Sol_HLOA,Convergence_curve_HLOA,NumberOfEval_HLOA] = HLOA(X,Fitness_Values,SearchAgents,Max_iterations,lowerbound,upperbound,dimension,fitness,AgentsReq,i);
        Best_Fitness_HLOA_nRuns(i,runs) = Best_Fitness_HLOA;
        Best_Sol_HLOA_nRuns(i,runs,:) = Best_Sol_HLOA;
        Convergence_curve_HLOA_nRuns(i,runs,:) = Convergence_curve_HLOA;
            
    end
end

disp([NumberOfEval_MSA_HO,NumberOfEval_MSA,NumberOfEval_HO,NumberOfEval_TLBO,NumberOfEval_HLOA]);

% Average Table
fprintf("\nAverage Table\n");
mean_case_1 = [abs(mean(Best_Fitness_MSA_HO_nRuns(1,:))),abs(mean(Best_Fitness_MSA_nRuns(1,:))),abs(mean(Best_Fitness_HO_nRuns(1,:))),abs(mean(Best_Fitness_TLBO_nRuns(1,:))),abs(mean(Best_Fitness_HLOA_nRuns(1,:)))];
mean_case_2 = [abs(mean(Best_Fitness_MSA_HO_nRuns(2,:))),abs(mean(Best_Fitness_MSA_nRuns(2,:))),abs(mean(Best_Fitness_HO_nRuns(2,:))),abs(mean(Best_Fitness_TLBO_nRuns(2,:))),abs(mean(Best_Fitness_HLOA_nRuns(2,:)))];
mean_case_3 = [abs(mean(Best_Fitness_MSA_HO_nRuns(3,:))),abs(mean(Best_Fitness_MSA_nRuns(3,:))),abs(mean(Best_Fitness_HO_nRuns(3,:))),abs(mean(Best_Fitness_TLBO_nRuns(3,:))),abs(mean(Best_Fitness_HLOA_nRuns(3,:)))];
mean_case_4 = [abs(mean(Best_Fitness_MSA_HO_nRuns(4,:))),abs(mean(Best_Fitness_MSA_nRuns(4,:))),abs(mean(Best_Fitness_HO_nRuns(4,:))),abs(mean(Best_Fitness_TLBO_nRuns(4,:))),abs(mean(Best_Fitness_HLOA_nRuns(4,:)))];
casess = ["case 1";"case 2";"case 3";"case 4"];
methods = ["Average_Table","Hybrid","MSA","HO","TLBO","HLOA"];
mean_result = [mean_case_1;mean_case_2;mean_case_3;mean_case_4];
mean_result = [casess,mean_result];
mean_result = [methods;mean_result];
disp(mean_result);

% Case 1 Stats Table
disp("Case 1 stats Table");
header1 = ["Case 1","Max","Min","Mean","Median","Std"];
col1 = ["Hybrid";"MSA";"HO";"TLBO";"HLOA"];
Hybrid_Stat = [abs(min(Best_Fitness_MSA_HO_nRuns(1,:))),abs(max(Best_Fitness_MSA_HO_nRuns(1,:))),abs(mean(Best_Fitness_MSA_HO_nRuns(1,:))),abs(median(Best_Fitness_MSA_HO_nRuns(1,:))),abs(std(Best_Fitness_MSA_HO_nRuns(1,:)))];
MSA_Stat    = [abs(min(Best_Fitness_MSA_nRuns(1,:))),abs(max(Best_Fitness_MSA_nRuns(1,:))),abs(mean(Best_Fitness_MSA_nRuns(1,:))),abs(median(Best_Fitness_MSA_nRuns(1,:))),abs(std(Best_Fitness_MSA_nRuns(1,:)))];
HO_Stat     = [abs(min(Best_Fitness_HO_nRuns(1,:))),abs(max(Best_Fitness_HO_nRuns(1,:))),abs(mean(Best_Fitness_HO_nRuns(1,:))),abs(median(Best_Fitness_HO_nRuns(1,:))),abs(std(Best_Fitness_HO_nRuns(1,:)))];
TLBO_Stat   = [abs(min(Best_Fitness_TLBO_nRuns(1,:))),abs(max(Best_Fitness_TLBO_nRuns(1,:))),abs(mean(Best_Fitness_TLBO_nRuns(1,:))),abs(median(Best_Fitness_TLBO_nRuns(1,:))),abs(std(Best_Fitness_TLBO_nRuns(1,:)))];
HLOA_Stat   = [abs(min(Best_Fitness_HLOA_nRuns(1,:))),abs(max(Best_Fitness_HLOA_nRuns(1,:))),abs(mean(Best_Fitness_HLOA_nRuns(1,:))),abs(median(Best_Fitness_HLOA_nRuns(1,:))),abs(std(Best_Fitness_HLOA_nRuns(1,:)))];

stat_1 = [Hybrid_Stat;MSA_Stat;HO_Stat;TLBO_Stat;HLOA_Stat];
stat_1 = [col1,stat_1];
stat_1 = [header1;stat_1];

disp(stat_1);

% Case 2 Stats Table
disp("Case 2 stats Table");
header2 = ["Case 2","Max","Min","Mean","Median","Std"];
col2 = ["Hybrid";"MSA";"HO";"TLBO";"HLOA"];
Hybrid_Stat = [abs(min(Best_Fitness_MSA_HO_nRuns(2,:))),abs(max(Best_Fitness_MSA_HO_nRuns(2,:))),abs(mean(Best_Fitness_MSA_HO_nRuns(2,:))),abs(median(Best_Fitness_MSA_HO_nRuns(2,:))),abs(std(Best_Fitness_MSA_HO_nRuns(2,:)))];
MSA_Stat    = [abs(min(Best_Fitness_MSA_nRuns(2,:))),abs(max(Best_Fitness_MSA_nRuns(2,:))),abs(mean(Best_Fitness_MSA_nRuns(2,:))),abs(median(Best_Fitness_MSA_nRuns(2,:))),abs(std(Best_Fitness_MSA_nRuns(2,:)))];
HO_Stat     = [abs(min(Best_Fitness_HO_nRuns(2,:))),abs(max(Best_Fitness_HO_nRuns(2,:))),abs(mean(Best_Fitness_HO_nRuns(2,:))),abs(median(Best_Fitness_HO_nRuns(2,:))),abs(std(Best_Fitness_HO_nRuns(2,:)))];
TLBO_Stat   = [abs(min(Best_Fitness_TLBO_nRuns(2,:))),abs(max(Best_Fitness_TLBO_nRuns(2,:))),abs(mean(Best_Fitness_TLBO_nRuns(2,:))),abs(median(Best_Fitness_TLBO_nRuns(2,:))),abs(std(Best_Fitness_TLBO_nRuns(2,:)))];
HLOA_Stat   = [abs(min(Best_Fitness_HLOA_nRuns(2,:))),abs(max(Best_Fitness_HLOA_nRuns(2,:))),abs(mean(Best_Fitness_HLOA_nRuns(2,:))),abs(median(Best_Fitness_HLOA_nRuns(2,:))),abs(std(Best_Fitness_HLOA_nRuns(2,:)))];

stat_2 = [Hybrid_Stat;MSA_Stat;HO_Stat;TLBO_Stat;HLOA_Stat];
stat_2 = [col2,stat_2];
stat_2 = [header2;stat_2];

disp(stat_2);

% Case 3 Stats Table
disp("Case 3 stats Table");
header3 = ["Case 3","Max","Min","Mean","Median","Std"];
col3    = ["Hybrid";"MSA";"HO";"TLBO";"HLOA"];
Hybrid_Stat = [abs(min(Best_Fitness_MSA_HO_nRuns(3,:))),abs(max(Best_Fitness_MSA_HO_nRuns(3,:))),abs(mean(Best_Fitness_MSA_HO_nRuns(3,:))),abs(median(Best_Fitness_MSA_HO_nRuns(3,:))),abs(std(Best_Fitness_MSA_HO_nRuns(3,:)))];
MSA_Stat    = [abs(min(Best_Fitness_MSA_nRuns(3,:))),abs(max(Best_Fitness_MSA_nRuns(3,:))),abs(mean(Best_Fitness_MSA_nRuns(3,:))),abs(median(Best_Fitness_MSA_nRuns(3,:))),abs(std(Best_Fitness_MSA_nRuns(3,:)))];
HO_Stat     = [abs(min(Best_Fitness_HO_nRuns(3,:))),abs(max(Best_Fitness_HO_nRuns(3,:))),abs(mean(Best_Fitness_HO_nRuns(3,:))),abs(median(Best_Fitness_HO_nRuns(3,:))),abs(std(Best_Fitness_HO_nRuns(3,:)))];
TLBO_Stat   = [abs(min(Best_Fitness_TLBO_nRuns(3,:))),abs(max(Best_Fitness_TLBO_nRuns(3,:))),abs(mean(Best_Fitness_TLBO_nRuns(3,:))),abs(median(Best_Fitness_TLBO_nRuns(3,:))),abs(std(Best_Fitness_TLBO_nRuns(3,:)))];
HLOA_Stat   = [abs(min(Best_Fitness_HLOA_nRuns(3,:))),abs(max(Best_Fitness_HLOA_nRuns(3,:))),abs(mean(Best_Fitness_HLOA_nRuns(3,:))),abs(median(Best_Fitness_HLOA_nRuns(3,:))),abs(std(Best_Fitness_HLOA_nRuns(3,:)))];

stat_3 = [Hybrid_Stat;MSA_Stat;HO_Stat;TLBO_Stat;HLOA_Stat];
stat_3 = [col3,stat_3];
stat_3 = [header3;stat_3];

disp(stat_3);

% Case 4 Stats Table
disp("Case 4 stats Table");

header4 = ["Case 4","Max","Min","Mean","Median","Std"];
col4 = ["Hybrid";"MSA";"HO";"TLBO";"HLOA"];

Hybrid_Stat = [abs(min(Best_Fitness_MSA_HO_nRuns(4,:))),abs(max(Best_Fitness_MSA_HO_nRuns(4,:))),abs(mean(Best_Fitness_MSA_HO_nRuns(4,:))),abs(median(Best_Fitness_MSA_HO_nRuns(4,:))),abs(std(Best_Fitness_MSA_HO_nRuns(4,:)))];
MSA_Stat    = [abs(min(Best_Fitness_MSA_nRuns(4,:))),abs(max(Best_Fitness_MSA_nRuns(4,:))),abs(mean(Best_Fitness_MSA_nRuns(4,:))),abs(median(Best_Fitness_MSA_nRuns(4,:))),abs(std(Best_Fitness_MSA_nRuns(4,:)))];
HO_Stat     = [abs(min(Best_Fitness_HO_nRuns(4,:))),abs(max(Best_Fitness_HO_nRuns(4,:))),abs(mean(Best_Fitness_HO_nRuns(4,:))),abs(median(Best_Fitness_HO_nRuns(4,:))),abs(std(Best_Fitness_HO_nRuns(4,:)))];
TLBO_Stat   = [abs(min(Best_Fitness_TLBO_nRuns(4,:))),abs(max(Best_Fitness_TLBO_nRuns(4,:))),abs(mean(Best_Fitness_TLBO_nRuns(4,:))),abs(median(Best_Fitness_TLBO_nRuns(4,:))),abs(std(Best_Fitness_TLBO_nRuns(4,:)))];
HLOA_Stat   = [abs(min(Best_Fitness_HLOA_nRuns(4,:))),abs(max(Best_Fitness_HLOA_nRuns(4,:))),abs(mean(Best_Fitness_HLOA_nRuns(4,:))),abs(median(Best_Fitness_HLOA_nRuns(4,:))),abs(std(Best_Fitness_HLOA_nRuns(4,:)))];

stat_4 = [Hybrid_Stat;MSA_Stat;HO_Stat;TLBO_Stat;HLOA_Stat];
stat_4 = [col4,stat_4];
stat_4 = [header4;stat_4];

disp(stat_4);

% Display best fitness and corresponding best solution for each process and case
for i = 1:cases
    fprintf("\nCase %d Best Solutions:\n", i);
    
    % Find the best fitness values and corresponding best solutions
    [best_fitness_Hybrid, best_idx_Hybrid] = min(Best_Fitness_MSA_HO_nRuns(i, :));
    best_solution_Hybrid = squeeze(Best_Sol_MSA_HO_nRuns(i, best_idx_Hybrid, :));
    
    [best_fitness_MSA, best_idx_MSA] = min(Best_Fitness_MSA_nRuns(i, :));
    best_solution_MSA = squeeze(Best_Sol_MSA_nRuns(i, best_idx_MSA, :));
    
    [best_fitness_HO, best_idx_HO] = min(Best_Fitness_HO_nRuns(i, :));
    best_solution_HO = squeeze(Best_Sol_HO_nRuns(i, best_idx_HO, :));
    
    [best_fitness_TLBO, best_idx_TLBO] = min(Best_Fitness_TLBO_nRuns(i, :));
    best_solution_TLBO = squeeze(Best_Sol_TLBO_nRuns(i, best_idx_TLBO, :));
    
    [best_fitness_HLOA, best_idx_HLOA] = min(Best_Fitness_HLOA_nRuns(i, :));
    best_solution_HLOA = squeeze(Best_Sol_HLOA_nRuns(i, best_idx_HLOA, :));
    
    % Display the best fitness and non-zero best solutions for each algorithm
    
    % Hybrid
    fprintf("\nHybrid: Best Fitness = %.2f, Best Solution (Non-zero) = [", abs(best_fitness_Hybrid));
    for j = 1:length(best_solution_Hybrid)
        % Check if value is above the 0.1 threshold to avoid very small values
        if abs(best_solution_Hybrid(j)) >= 0.1
            fprintf("\nProcess %d: %.2f , Product - %d", j, best_solution_Hybrid(j),product(j));
        end
    end
    fprintf("]\n");
    
    % MSA
    fprintf("\nMSA: Best Fitness = %.2f, Best Solution (Non-zero) = [", abs(best_fitness_MSA));
    for j = 1:length(best_solution_MSA)
        if abs(best_solution_MSA(j)) >= 0.1
            fprintf("\nProcess %d: %.2f, Product - %d ", j, best_solution_MSA(j),product(j));
        end
    end
    fprintf("]\n");
    
    % HO
    fprintf("\nHO: Best Fitness = %.2f, Best Solution (Non-zero) = [", abs(best_fitness_HO));
    for j = 1:length(best_solution_HO)
        if abs(best_solution_HO(j)) >= 0.1
            fprintf("\nProcess %d: %.2f, Product - %d  ", j, best_solution_HO(j),product(j));
        end
    end
    fprintf("]\n");
    
    % TLBO
    fprintf("\nTLBO: Best Fitness = %.2f, Best Solution (Non-zero) = [", abs(best_fitness_TLBO));
    for j = 1:length(best_solution_TLBO)
        if abs(best_solution_TLBO(j)) >= 0.1
            fprintf("\nProcess %d: %.2f, Product - %d  ", j, best_solution_TLBO(j),product(j));
        end
    end
    fprintf("]\n");
    
    % HLOA
    fprintf("\nHLOA: Best Fitness = %.2f, Best Solution (Non-zero) = [", abs(best_fitness_HLOA));
    for j = 1:length(best_solution_HLOA)
        if abs(best_solution_HLOA(j)) >= 0.1
            fprintf("\nProcess %d: %.2f, Product - %d  ", j, best_solution_HLOA(j),product(j));
        end
    end
    fprintf("]\n");
end



% Plot convergence curves for each case
figure;
for i = 1:cases
    subplot(2, 2, i); % Set up a 2x2 grid of subplots
    
    % Retrieve convergence curves for each method for the current case
    curve_Hybrid = squeeze(Convergence_curve_MSA_HO_nRuns(i, :, :));
    curve_MSA = squeeze(Convergence_curve_MSA_nRuns(i, :, :));
    curve_HO = squeeze(Convergence_curve_HO_nRuns(i, :, :));
    curve_TLBO = squeeze(Convergence_curve_TLBO_nRuns(i, :, :));
    curve_HLOA = squeeze(Convergence_curve_HLOA_nRuns(i, :, :));
    
    % Find the best convergence curve (minimum final fitness value) for each method
    [~, best_idx_Hybrid] = min(Best_Fitness_MSA_HO_nRuns(i, :));
    [~, best_idx_MSA] = min(Best_Fitness_MSA_nRuns(i, :));
    [~, best_idx_HO] = min(Best_Fitness_HO_nRuns(i, :));
    [~, best_idx_TLBO] = min(Best_Fitness_TLBO_nRuns(i, :));
    [~, best_idx_HLOA] = min(Best_Fitness_HLOA_nRuns(i, :));
    
    % Inline code to find the starting index for negative values for each method
    start_Hybrid = find(curve_Hybrid(best_idx_Hybrid, :) < 0, 1);
    if isempty(start_Hybrid), start_Hybrid = length(curve_Hybrid(best_idx_Hybrid, :)); end
    
    start_MSA = find(curve_MSA(best_idx_MSA, :) < 0, 1);
    if isempty(start_MSA), start_MSA = length(curve_MSA(best_idx_MSA, :)); end
    
    start_HO = find(curve_HO(best_idx_HO, :) < 0, 1);
    if isempty(start_HO), start_HO = length(curve_HO(best_idx_HO, :)); end
    
    start_TLBO = find(curve_TLBO(best_idx_TLBO, :) < 0, 1);
    if isempty(start_TLBO), start_TLBO = length(curve_TLBO(best_idx_TLBO, :)); end
    
    start_HLOA = find(curve_HLOA(best_idx_HLOA, :) < 0, 1);
    if isempty(start_HLOA), start_HLOA = length(curve_HLOA(best_idx_HLOA, :)); end
    
    % Plot convergence curves for the best run of each method from the first negative value onwards
    plot(curve_Hybrid(best_idx_Hybrid, start_Hybrid:end), 'DisplayName', 'Hybrid', 'LineWidth', 1.5); hold on;
    plot(curve_MSA(best_idx_MSA, start_MSA:end), 'DisplayName', 'MSA', 'LineWidth', 1.5);
    plot(curve_HO(best_idx_HO, start_HO:end), 'DisplayName', 'HO', 'LineWidth', 1.5);
    plot(curve_TLBO(best_idx_TLBO, start_TLBO:end), 'DisplayName', 'TLBO', 'LineWidth', 1.5);
    plot(curve_HLOA(best_idx_HLOA, start_HLOA:end), 'DisplayName', 'HLOA', 'LineWidth', 1.5);
    
    % Labels and title for each subplot
    xlabel('Iterations');
    ylabel('Fitness Value');
    title(['Convergence Plot - Case ', num2str(i)]);
    legend('show'); % Show legend for each subplot
    hold off;
end

