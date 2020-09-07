%Use to compute the COMIC and AIC of the random walk
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
intNumOfEns = 30; %number of ensemble for random walk
Domain = [-5 5] + dblXShift;
dblBinSize = 0.1;

vecNumData = [10 30 200];
vecColor = {'r','b','k'};
vecdx = 10./(vecNumData - 1);

matExactSolution = zeros(length(vecNumData),vecNumData(end));
matX = zeros(length(vecNumData),vecNumData(end));
for i = 1:length(vecNumData)
    matX(i,1:vecNumData(i)) = Domain(1,1):vecdx(i):Domain(1,2);
    %matX(i,1:vecNumData(i)) = sort(10*rand(1,vecNumData(i)) - 5);
    matExactSolution(i,1:vecNumData(i)) = ExactSolution1D(matX(i,1:vecNumData(i)),T,D,V,xCenter);
end


intIter = 12;
vecIter = 1:1:intIter;
vecNumOfParts = 25*2.^vecIter;

matAICResults = zeros(intIter,3);
celLegend = cell(1,length(vecNumData)*2);
intNumOfTrial = 200;
DValue = D;

figure(1)
for k = 1:length(vecNumData)
    
    for j=1:intNumOfTrial
        
        funSSEN =@(n)1/vecNumData(k)*norm(matExactSolution(k,1:vecNumData(k))'-...
            ApproxSolution1D(matX(k,1:vecNumData(k)),dblBinSize,Domain,DValue,T,n,intNumOfEns,V,xCenter,'B'),2)^2;
        
        for i=1:intIter
            matAICResults(i,1) = vecNumOfParts(i);
%            if vecNumData(k) < 20
                matAICResults(i,2) = (matAICResults(i,2)*(j-1) + 2*log(funSSEN(vecNumOfParts(i))) + 12/(vecNumData(k)- 3))/j;
%            else
%                 matAICResults(i,2) = (matAICResults(i,2)*(j-1) + 2*log(funSSEN(vecNumOfParts(i))))/j;
%            end
%            matAICResults(i,2) = (matAICResults(i,2)*(j-1) + 2*log(funSSEN(vecNumOfParts(i))))/j;
            matAICResults(i,3) = (matAICResults(i,3)*(j-1) + matAICResults(i,2) + log(vecNumOfParts(i)))/j;
        end
        
    end
    plot(log10(matAICResults(:,1)),matAICResults(:,2),'-x','color',string(vecColor(k)));
    hold on
    plot(log10(matAICResults(:,1)),matAICResults(:,3),'-o','color',string(vecColor(k)))
    celLegend{1,k*2-1} = strcat('AIC-',num2str(vecNumData(k)));
    celLegend{1,k*2} = strcat('COMIC-',num2str(vecNumData(k)));
    hold on
    k
end
legend(celLegend,'Location','southwest','Interpreter', 'latex')
xlabel('$\log_{10}$($n$)','Interpreter', 'latex')
ylabel('Fitness Metric','Interpreter', 'latex')
hold off

