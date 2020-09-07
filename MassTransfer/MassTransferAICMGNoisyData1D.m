%Use to compute the COMIC and AIC with non-Gaussian error process where the
%data actually have some noise
clc; clear all;
close all;
rng('shuffle')
set(0,'defaultTextInterpreter','latex');
%Initialize
D = 1;
T = 1;
V = 0;
xCenter = 0;
dblXShift = V*T + xCenter;
intNumData = 30;%number of data points
intNumOfPart = 3000; %number of particles
intNumOfEns = 1; %number of ensemble for random walk
Domain = [-6 6];
ZeroDomain = Domain;
Domain = Domain + dblXShift;
%dx = (Domain(1,2) - Domain(1,1))/intNumData;
dx = 10/(intNumData - 1);
dblBinSize = 0.1;
%dblAlpha = 1/10;
dblAlpha = 1;
%Choose x points
%x = (Domain(1,2) - Domain(1,1))*rand(1,intNumData) + Domain(1,1);
%x = 10*rand(1,intNumData) - 5;
%x = Domain(1,1):dx:Domain(1,2);
x = -5:dx:5;
x = sort(x);
intNumData = length(x);
vecErrorNoise = [3 9 81];
celLegend = cell(1,length(vecErrorNoise*2));
celColor = {'r','b','k'};

for k=1:length(vecErrorNoise)
    %Get exact solution
    vecExactSolution = ExactSolution1D(x,T,D,V,xCenter);
    vecTest = vecExactSolution;
    %Set up the error to the data
    dblDelta = min(vecExactSolution);
    for s=1:length(vecExactSolution)
        distNNoise = makedist('Normal');
        distNNoise.mu = 0;
        distNNoise.sigma = vecExactSolution(s)/vecErrorNoise(k);
        vecExactSolution(s) = vecExactSolution(s) + random(distNNoise);
    end
    
    %Find D that minimize the error
    % funErrorD =@(d) 1/intNumData*norm(vecExactSolution'-...
    %     MassTransferImp1DSolution(x,d,Domain,ZeroDomain,intNumOfPart,T,intNumOfEns,dx,V,xCenter),2)^2;
    funErrorD =@(d) 1/intNumData*sum(1./(dblAlpha.*vecExactSolution').*(vecExactSolution'-...
        MassTransferImp1DSolution(x,d,Domain,ZeroDomain,intNumOfPart,T,intNumOfEns,dblBinSize,V,xCenter)).^2);
    DValue = 0.5;%inital guess
    
    %----------------------------------------------------------------------
    %These 2 line are alternative of each other, allow to use the estimate
    %or the actual value, modify to see how thing behave. Need to comment
    %one out for the other one to work
    %[DValue, DError] = fminsearch(funErrorD,DValue);
    DValue = 1;
    
    %---------------------------------------------------------------------
    %Using this approx D value to compute the AIC and AIC + ln(N)
    vecNumOfParts = [20 30 100 300 500 1000 3000 8000 18000]';
    intIter = length(vecNumOfParts);
    matAICResults = zeros(intIter,3);
    
    intNumOfTrial = 1;

    dblAlpha = 1;
    for j=1:intNumOfTrial
        
        for i=1:intIter
            matAICResults(i,1) = vecNumOfParts(i);
            dblError = 1/intNumData*sum(1./(dblAlpha.*vecExactSolution').*(vecExactSolution'-...
                MassTransferImp1DSolution(x,DValue,Domain,ZeroDomain,vecNumOfParts(i),T,intNumOfEns,dblBinSize,V,xCenter)).^2);
            matAICResults(i,2) = (matAICResults(i,2)*(j-1) + 2*log(dblError))/j;
            
            matAICResults(i,3) = (matAICResults(i,3)*(j-1) + matAICResults(i,2) + log(vecNumOfParts(i)))/j;
        end
        
    end
    plot(log10(matAICResults(:,1)),matAICResults(:,2),'-x','color',string(celColor(k)));
    hold on
    plot(log10(matAICResults(:,1)),matAICResults(:,3),'-o','color',string(celColor(k)))
    
    celLegend{1,k*2-1} = strcat('$2\ln(\mathcal{E}):\alpha =1/$',num2str(vecErrorNoise(k)));
    celLegend{1,k*2} = strcat('COMIC:$\alpha=1/$',num2str(vecErrorNoise(k)));
    hold on
end

legend(celLegend,'Location','southwest','Interpreter', 'latex')
xlabel('$\log_{10}(n)$','Interpreter', 'latex')
ylabel('Fitness Metric','Interpreter', 'latex')
ylim([-35 5]);


