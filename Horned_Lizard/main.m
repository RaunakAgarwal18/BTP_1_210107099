%  ********************************************
%  HLOA: Horned Lizard Optimization Algorithm 
%  Developed in MATLAB R2018a(9.4)                                        %
%  Author/Programmer:  Dr. Hernan Peraza-Vazquez                                     %
%  email: hperaza@ipn.mx                                %
%  Telegram: @CodebugMx                                                                   %
%  Paper: A Novel Metaheuristic Inspired by Horned Lizard Defense Tactics %
%  Artificial Intelligence Review - Springer (2024) 
%  Doi: 10.1007/s10462-023-10653-7
%_________________________________________________________________________%
clear
clc
format long
i=5;  % number of executions
a=zeros(i,1);
SearchAgents_no=30;  
Function_name= 0;
Max_iteration=50; 
[lb,ub,dim,fobj]=GetFunction(Function_name);
v= zeros(i,dim);
cc= zeros(i,Max_iteration);   % convergence curve
for m=1:i
   % display(['Run: ', num2str(m)]);
   [vMin,theBestVct,Convergence_curve]=HLOA(SearchAgents_no,Max_iteration,lb,ub,dim,fobj); 
   a(m)=vMin;  % the best fitness for the m-Run 
   v(m,:)= theBestVct; % the best vector (variables) for the m-Run
   cc(m,:)= Convergence_curve; % convergence curve for the m-Run
end
[vMin,minIdx]= min(a); 
theBestVct = v(minIdx,:);
ConvergenceC= cc(minIdx,:);

disp(min(ConvergenceC));
% Plot the convergence curve for the best run
figure;
plot(ConvergenceC, 'LineWidth', 2);
axis auto;
title('Convergence Curve');
xlabel('Iteration');
ylabel('Fitness Value');
% grid on;
% axis([0 800 -600 1e24]);


display(['The best solution obtained by HLOA is: ', num2str(theBestVct)]);
display(['The best fitness (min f(x)) found by HLOA is: ', num2str(vMin)]);
display(['# runs: ', num2str(i)]);
display(['Mean:  ', num2str(mean(a))]);
display(['Std.Dev:  ', num2str(std(a))]);
%%*************************************************

