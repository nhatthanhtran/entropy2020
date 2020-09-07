%Use to calculate the fitness metric using the estimate parameters at each
%number of particles.

clc; clear all; 
close all;
rng('shuffle')
set(0,'defaultTextInterpreter','latex');
%Initialize
D = 1;
T = 1;
V = 1;
xCenter = 0;
dblXShift = V*T + xCenter;
intNumData = 30;%number of data points
intNumOfPart = 300; %number of particles
intNumOfEns = 1; %number of ensemble for random walk
Domain = [-6 6];
ZeroDomain = Domain;
Domain = Domain + dblXShift;
%dx = (Domain(1,2) - Domain(1,1))/intNumData; 
dx = 10/(intNumData - 1);
dblBinSize = 0.1;


x = -4:dx:6;

x = sort(x);
intNumData = length(x);
%Get exact solution
vecExactSolution = ExactSolution1D(x,T,D,V,xCenter);

%vecNumOfParts = [28 30 32 34 36 38 40 100 300 500 1000 2000 3000 4000 5000 6000]';
vecNumOfParts = [30 100 300 500 1000 3000 5000 8000]';
intIter = length(vecNumOfParts);
matAICResults = zeros(intIter,3);
matDVValue = zeros(intIter,2);
intNumOfTrial = 1;
%DValue = 1;
matError = zeros(intNumData,intIter);

tic
for j=1:intNumOfTrial
    
    for i=1:intIter
        funErrorD =@(dv) 1/intNumData*norm(vecExactSolution'-...
            MassTransferImp1DSolution(x,dv(1),Domain,ZeroDomain,vecNumOfParts(i),T,intNumOfEns,dx,dv(2),xCenter),2)^2;
        DValue = [0.5 0.5];%inital guess
        
        %Estimate the parameters
        [DValue, DError] = fminsearch(funErrorD,DValue);
        
        
        matDVValue(i,:) = DValue;
        matAICResults(i,1) = vecNumOfParts(i);
        dblError = 1/intNumData*norm(vecExactSolution'-...
            MassTransferImp1DSolution(x,DValue(1),Domain,ZeroDomain,vecNumOfParts(i),T,intNumOfEns,dblBinSize,DValue(2),xCenter),2)^2;
        matAICResults(i,2) = (matAICResults(i,2)*(j-1) + 2*log(dblError) + 12/(intNumData - 3))/j;
        matAICResults(i,3) = (matAICResults(i,3)*(j-1) + matAICResults(i,2) + log(vecNumOfParts(i)))/j;
        i
    end
    
end
figure(1)
plot(log10(matAICResults(:,1)),matAICResults(:,2),'-x',log10(matAICResults(:,1)),matAICResults(:,3),'-o');
legend('AIC', 'COMIC','Interpreter', 'latex')
xlabel('$\log_{10}$($n$)','Interpreter', 'latex')
ylabel('Fitness Metric','Interpreter', 'latex')
toc
