%____________________ 1. Reconnaissance phase : Pathfinder moths ____________________________________________________    
function [Lights,Light_fitness] = Pathfinders(Lights,Light_fitness,Nc,fobj,ub,lb,caseNo)

Lights1 = Lights;  %trail vector

%Proposed adaptive crossover operation based on population diversity
% coeficient of variation
C_V  = std(Lights,0,1)./abs(mean(Lights));             % sensing_distance using Eqs.(18) and (19)
nmu  = mean(C_V);
mcol = find(C_V <= nmu);                               % Select crossover points using condition of Eq.(20)

%levy mutations for crossover using Eq.(24)
beta  = 3/2;
sigma = (gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u     = randn(Nc,length(mcol),2)*sigma;
v     = randn(Nc,length(mcol),2);
step  = u./abs(v).^(1/beta);
L     = abs(0.01*step);                                % Draw Levy flight sample

% Proposed difference vectors Lévy-mutation 
for i = 1:Nc
    if Nc<6
        I               = randperm(Nc);  
        II              = find(I ~= i);
        I               = I(II(1:3));
        Lights1(i,mcol) = Lights(I(1),mcol)+L(i,:,1).*(Lights(I(2),mcol)-Lights(I(3),mcol));
    else
        I               = randperm(Nc);  
        II              = find(I ~= i);
        I               = I(II(1:5));   %mating vectors (donners)
        Lights1(i,mcol) = Lights(I(1),mcol)+(L(i,:,1).*(Lights(I(2),mcol)-Lights(I(3),mcol))+L(i,:,2).*(Lights(I(4),mcol)-Lights(I(5),mcol))); 
    end
end

% Return back the Pathfinders that go beyond the boundaries
Lights1 = Bound_Checking(Lights1,ub,lb);

%selection for updating Pathfinder moths using Eq.(27)
for i = 1:size(Lights1,1)
    [Lights1(i,:),Light_fitness1(i,:)] = fobj(Lights1(i,:),caseNo); %#ok<AGROW>
    if Light_fitness1(i,1) < Light_fitness(i,1)
        Light_fitness(i,:) = Light_fitness1(i,:);
        Lights(i,:)        = Lights1(i,:);
    end
end

