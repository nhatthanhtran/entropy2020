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
intNumData = 30; %number of data points
intNumOfPart = 5000; %number of particles
intNumOfEns = 15; %number of ensemble for random walk
Domain = [-5 5] + dblXShift;
dx = (Domain(1,2) - Domain(1,1))/(intNumData-1); 
dblBinSize = 0.1;
%Choose x points
%x = (Domain(1,2) - Domain(1,1))*rand(1,intNumData) + Domain(1,1);
x = Domain(1,1):dx:Domain(1,2);
x = sort(x);
intNumData = length(x);
vecErrorNoise = [3 9 81];
celLegend = cell(1,length(vecErrorNoise*2));
celColor = {'r','b','k'};
for k = 1:length(vecErrorNoise)
%Get exact solution
vecExactSolution = ExactSolution1D(x,T,D,V,xCenter);
vecTest  =  vecExactSolution;
vecMu = zeros(length(vecTest),2);
for s=1:length(vecExactSolution)
    distNNoise = makedist('Normal');
    distNNoise.mu = 0;
    distNNoise.sigma = vecExactSolution(s)/vecErrorNoise(k);
    vecMu(s,1) = vecExactSolution(s);
    vecMu(s,2) = distNNoise.sigma;

    vecExactSolution(s) = vecExactSolution(s) + random(distNNoise); 
end
if min(vecExactSolution) < 0
   disp("Bad Data") 
end
%Find D that minimize the error
% funErrorD =@(d) 1/intNumData*norm(vecExactSolution'-...
%     ApproxSolution1D(x,dblBinSize,Domain,d,T,intNumOfPart,intNumOfEns,V,xCenter,'B'),2)^2;

funErrorD = @(d)1/intNumData*sum(1./vecExactSolution'.*(vecExactSolution'- ...
            ApproxSolution1D(x,dblBinSize,Domain,d,T,intNumOfPart,intNumOfEns,V,xCenter,'B')).^2);

DValue = 0.5;%inital guess
%--------------------------------------------------------------------------
%Use to set up the D value in which there maybe an estimate or actual
%[DValue, DError] = fminsearch(funErrorD,DValue);
DValue = 1;
%--------------------------------------------------------------------------
%Using this approx D value to compute the AIC and AIC + ln(N)
intIter = 12;
vecIter = 1:1:intIter;
vecNumOfParts = [5*2.^vecIter];

matAICResults = zeros(intIter,3);

intNumOfTrial = 30;

for j=1:intNumOfTrial
    
    for i=1:intIter
        dblAlpha = 1;
        funSSEN =@(n) 1/intNumData*sum(1./(dblAlpha.*vecExactSolution').*(vecExactSolution'- ...
            ApproxSolution1D(x,dblBinSize,Domain,D,T,n,intNumOfEns,V,xCenter,'B')).^2);
        matAICResults(i,1) = vecNumOfParts(i);
        matAICResults(i,2) = (matAICResults(i,2)*(j-1) + 2*log(funSSEN(vecNumOfParts(i))))/j;
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
hold off
