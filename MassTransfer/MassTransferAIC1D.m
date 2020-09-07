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
intNumOfEns = 1; %number of ensemble if needed, no longer need
Domain = [-6 6];
ZeroDomain = Domain;
Domain = Domain + dblXShift;
dblBinSize = 0.1;%no longer need

%Set up the number of data
vecNumData = [10 30 200];
celColor = {'r','b','k'};
vecdx = 10./(vecNumData - 1);

%Set up all of the exact solutions
matExactSolution = zeros(length(vecNumData),vecNumData(end));
matX = zeros(length(vecNumData),vecNumData(end));
for i = 1:length(vecNumData)
    matX(i,1:vecNumData(i)) = -5:vecdx(i):5;
    %matX(i,1:vecNumData(i)) = sort(10*rand(1,vecNumData(i)) - 5);
    matExactSolution(i,1:vecNumData(i)) = ExactSolution1D(matX(i,1:vecNumData(i)),T,D,V,xCenter);
end

%Set up the number of particles
vecNumOfParts = [30 100 300 500 1000 3000 5000 8000 18000];
intIter = length(vecNumOfParts);

matAICResults = zeros(intIter,3);

intNumOfTrial = 1;
DValue = D;
celLegend = cell(1,length(vecNumData)*2);
clf
figure(1)
for k = 1:length(vecNumData)
    
    for j=1:intNumOfTrial
        for i=1:intIter
            matAICResults(i,1) = vecNumOfParts(i);
            dblError = 1/vecNumData(k)*norm(matExactSolution(k,1:vecNumData(k))'-...
                MassTransferImp1DSolution(matX(k,1:vecNumData(k)),DValue,Domain,ZeroDomain,vecNumOfParts(i),T,intNumOfEns,dblBinSize,V,xCenter),2)^2;
            matAICResults(i,2) = (matAICResults(i,2)*(j-1) + 2*log(dblError) + 12/(vecNumData(k)- 3))/j;
            matAICResults(i,3) = (matAICResults(i,3)*(j-1) + matAICResults(i,2) + log(vecNumOfParts(i)))/j;
            
        end
        
    end
    plot(log10(matAICResults(:,1)),matAICResults(:,2),'-x','color',string(celColor(k)));
    hold on
    plot(log10(matAICResults(:,1)),matAICResults(:,3),'-o','color',string(celColor(k)))
    celLegend{1,k*2-1} = strcat('AIC-',num2str(vecNumData(k)));
    celLegend{1,k*2} = strcat('COMIC-',num2str(vecNumData(k)));
    
    hold on
end
legend(celLegend,'Location','southwest','Interpreter', 'latex')
xlabel('$\log_{10}$($n$)','Interpreter', 'latex')
ylabel('Fitness Metric','Interpreter', 'latex')
hold off



