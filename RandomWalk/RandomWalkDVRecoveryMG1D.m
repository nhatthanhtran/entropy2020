%Parameter estimation for random walk - using the non Gaussian process
clc; clear all; 
close all;
rng('shuffle');
set(0,'defaultTextInterpreter','latex');
%Try to recover D from a noisy data
%Initialize
D = 1;
T = 1;
V = 1;
xCenter = 0;
dblXShift = V*T + xCenter;
intNumData = 30;%number of data points
intNumOfPartStar = 5000;%NEED to update this accordingly to the number of dataa
intNumOfEns = 15; %number of ensemble for random walk
Domain = [-5 5] + dblXShift;
dx = (Domain(1,2) - Domain(1,1))/(intNumData-1); 
dblBinSize = 0.1;

%Randomly choosing the x data point within the domain
%x = (Domain(1,2) - Domain(1,1))*rand(1,intNumData) + Domain(1,1);
x = Domain(1,1):dx:Domain(1,2);
x = sort(x);
intNumData = length(x);

%Get exact solution
vecExactSolution = ExactSolution1D(x,T,D,V,xCenter);

dblAlpha = 1;

intNumOfTrials = 30;
matResults = zeros(intNumOfTrials,2);
for i =1:intNumOfTrials
   funErrorDStar =@(d) 1/intNumData*sum(1./(dblAlpha.*vecExactSolution').*(vecExactSolution'- ...
        ApproxSolution1D(x,dblBinSize,Domain,d(1),T,intNumOfPartStar,intNumOfEns,d(2),xCenter,'B')).^2);
   DValueStar = [0.5;0.5]; 
   [DValueStar, DErrorStar] = fminsearch(funErrorDStar,DValueStar);
   matResults(i,:) = DValueStar';
end

figure(1)
subplot(1,2,1)
boxplot(matResults(:,1),'Labels',{''})
ylabel('Estimate $D$', 'Interpreter', 'latex')
subplot(1,2,2)
boxplot(matResults(:,2),'Labels',{''})
ylabel('Estimate $v$', 'Interpreter', 'latex')
