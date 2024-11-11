function [x,f] = SKS_ProductionPlanningC(x,caseNo)

[product,l,m,h,il,im,ih,cl,cm,ch,SalePrice,rm1,rm2,~] = ProductionPlanningData; % Accessing the data from the function
nProcess = length(l);                         % determing the number of processes

if caseNo == 1
    Budget    = 1000;
    AvailRaw1 = 500;
    AvailRaw2 = 500;
elseif caseNo == 2
    Budget    = 1000;
    AvailRaw1 = 1000;
    AvailRaw2 = 1000;
elseif caseNo == 3
    Budget    = 2000;
    AvailRaw1 = 500;
    AvailRaw2 = 500;
else
    Budget    = 2000;
    AvailRaw1 = 1000;
    AvailRaw2 = 1000;
end

PC = zeros(nProcess,1);                     % Initialization of variable to store the production cost for each process
IC = zeros(nProcess,1);                     % Initialization of variable to store the investment cost for each process
R1Reqd = zeros(nProcess,1);                 % Initialization of variable to store the raw material 1 required for each process
R2Reqd = zeros(nProcess,1);                 % Initialization of variable to store the raw material 2 required for each process
% penalty_domain = zeros(nProcess,1);       % Initialization of variable to store the penalty due to violation of domain for each process
Revenue = zeros(nProcess,1);                % Initialization of variable to store the revenue genereated for each process




% Unique Process Constraints Implementation
startIdx = 1;
endIdx = 1;
index = 1;
indexes = [];  % Initialize indexes array to store active process indices

while endIdx <= nProcess
    if(product(endIdx) == product(startIdx))
        % Check if x(endIdx) lies within the limits l(endIdx) and h(endIdx)
        if x(endIdx) >= l(endIdx) && x(endIdx) <= h(endIdx)
            indexes(index) = endIdx;  % Store the index of active processes
            index = index + 1;
        end
        endIdx = endIdx + 1;  % Move to the next process
    else
        % Choose a random index from active processes in the current product group
        if ~isempty(indexes)
            random_index = indexes(randi(length(indexes)));  % Select a random element from indexes
        end

        % Set all other processes for this product to 0
        for i = indexes
            if(i ~= random_index)
                x(i) = 0;
            end
        end

        % Reset for the next group
        startIdx = endIdx;
        indexes = [];  % Reset the active process index list
        index = 1;
    end
end

% Final check in case last group reached end without exiting the loop
if ~isempty(indexes)
    random_index = indexes(randi(length(indexes)));
    for i = indexes
        if(i ~= random_index)
            x(i) = 0;
        end
    end
end





%TLBO Starts

for j = 1: nProcess
    
    if x(j) >= l(j) && x(j) <= m(j)             % if the production is in between low and medium level
        
        PC(j) = ((cm(j) - cl(j))/(m(j) - l(j)))*(x(j) - l(j)) + cl(j);  % Determining the production cost for the range l to m
        IC(j) = ((im(j) - il(j))/(m(j) - l(j)))*(x(j) - l(j)) + il(j);  % Determining the investment cost  for the range l to m
        R1Reqd(j) = x(j)*rm1(j);                % Determining the raw material 1 required for process
        R2Reqd(j) = x(j)*rm2(j);                % Determining the raw material 2 required for process
        Revenue(j) = SalePrice(j)*x(j);         % Determining the revenue genereated by selling the product from the process j
        
        
    elseif x(j) > m(j) && x(j) <= h(j)          % if the production is in between medium and high level
        
        PC(j) = ((ch(j) - cm(j))/(h(j) - m(j)))*(x(j) - m(j)) + cm(j);  % Determining the production cost
        IC(j) = ((ih(j) - im(j))/(h(j) - m(j)))*(x(j) - m(j)) + im(j);  % Determining the investment cost
        R1Reqd(j) = x(j)*rm1(j);                 % Determining the raw material 1 required for process
        R2Reqd(j) = x(j)*rm2(j);                 % Determining the raw material 2 required for process
        Revenue(j) = SalePrice(j)*x(j);          % Determining the revenue genereated by selling the product from the process j
        
        
    elseif 0 < x(j) &&  x(j)< l(j)                % if the production is greater than zero but less than the low level
        x(j) = 0;
        % penalty_domain(j) = 10^5;               % Assigning penalty for violating the domain hole constraints
    end
    
end

TotalInvRed = sum(IC);                           % Total Investment Cost Required
penalty_IC = 0;
if  TotalInvRed > Budget                         % Checking for violation of the investment cost constraint
    penalty_IC = (TotalInvRed - Budget)^2;       % Determining the penalty due to violation for investment cost
end

penalty_R1 = 0;
TotalR1Reqd = sum(R1Reqd);
if TotalR1Reqd > AvailRaw1                       % Checking for violation of the raw material 1 constraint
    penalty_R1 = (TotalR1Reqd  - AvailRaw1)^2;   % Determining the penalty due to violation for raw material 1
end

penalty_R2 = 0;
TotalR2Reqd = sum(R2Reqd);
if TotalR2Reqd > AvailRaw2                       % Checking for violation of the raw material 2 constraint
    penalty_R2 = (TotalR2Reqd  - AvailRaw2)^2;   % Determining the penalty due to violation for raw material 2
end

profit = sum(Revenue) - sum(PC);                 % Determining the profit (objective)

f = -profit + 10^15*(penalty_IC + penalty_R1 + penalty_R2);   % Determining the fitness function value

% + sum(penalty_domain)