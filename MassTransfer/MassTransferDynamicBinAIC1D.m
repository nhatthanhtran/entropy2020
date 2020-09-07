%Use to compute the COMIC and AIC of the non uniform binning for mass
%transfer
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
intNumOfPart = 300; %number of particles
intNumOfEns = 1; %number of ensemble for random walk
Domain = [-6 6];
ZeroDomain = Domain;
Domain = Domain + dblXShift;
%dx = (Domain(1,2) - Domain(1,1))/intNumData;
dx = (Domain(1,2)-Domain(1,1))/(intNumData - 1);
dblBinSize = 0.1;
%Choose x points
%x = (Domain(1,2) - Domain(1,1))*rand(1,intNumData) + Domain(1,1);
%x = Domain(1,1):dx:Domain(1,2);
x = -5:10/(intNumData-1):5;
x = sort(x);
intNumData = length(x);
%Get exact solution
vecExactSolution = ExactSolution1D(x,T,D,V,xCenter);

vecNumOfParts = [30 100 300 500 1000 3000 5000 8000 18000 25000 30000]';
%vecNumOfParts = [10000 40000]
intIter = length(vecNumOfParts);
matAICResults = zeros(intIter,4);
vecFineGrid = Domain(1,1):0.01:Domain(1,2);%This is use to calculate the actual 
vecFineGrid = -5:0.01:5;
%integral of the COMIC

intNumOfTrial = 30;
DValue = 1;%Set up the D Value, this can modify to see the effect and also 
%if one want to to parameter estimation before computation
tic
for j=1:intNumOfTrial
    
    for i=1:intIter
        %Set up the initial particle position, this need to be done to
        %avoid alot of spiking in the solution, this remove irregularity
        %to the solution -smooth out the solution to have better fit
        matAICResults(i,1) = vecNumOfParts(i);
        intNumOfParts = vecNumOfParts(i);
        dblNormalSpace = 12/(intNumOfParts-1);
        dblCenter = V*T + xCenter;
        dblLeftPart = dblCenter - dblNormalSpace/2;
        dblRightPart = dblCenter + dblNormalSpace/2;
        dblLLeftPart = dblLeftPart - dblNormalSpace;
        dblRRightPart = dblRightPart + dblNormalSpace;
        
        intLNumOfParts = (intNumOfParts - 4)/2;
        intRNumOfParts = (intNumOfParts - 4)/2;
        vecLeftBin = [Domain(1,1) (dblLLeftPart - Domain(1,1))*rand(1,intLNumOfParts-1) + Domain(1,1)];
        vecRightBin = [(Domain(1,2) - dblRRightPart)*rand(1,intRNumOfParts-1) + dblRRightPart Domain(1,2)];
        vecBin =[vecLeftBin dblLLeftPart dblLeftPart dblRightPart dblRRightPart vecRightBin];
        
        clear vecLeftBin vecRightBin
        %vecBin = Domain(1,1):((Domain(1,2)-Domain(1,1))/(intNumOfParts-1)):Domain(1,2);
        vecBin = sort(vecBin);
        
        [vecApproxSolution,matBin,vecCon,intCenterIdx] = MassTransferImpDynamicBin1DSolution(x,DValue,Domain,ZeroDomain,intNumOfParts,T,intNumOfEns,vecBin,V,xCenter);
        
        dblInformation = 0;
        dblAIC = 1/intNumData*norm(vecExactSolution'-vecApproxSolution,2)^2;
        
        %Calculating the actual COMIC with the integral formula
        p = 1;
        u = 1;
        i

        for k = 1:length(vecFineGrid)-1
            %Find the bin
            for q=p:length(vecBin)-1
                if vecFineGrid(k) >= vecBin(q) && vecFineGrid(k) <=vecBin(q+1)
                    p = q;
                    break;
                end
            end
            
            failcount = 0;
            if vecCon(p) > 1e-8
                dblBinSize = vecBin(p+1) - vecBin(p);
                
                dbldx = vecFineGrid(k+1) - vecFineGrid(k);
                dblIValue = vecCon(p)*log(dblBinSize)*dbldx;
                dblInformation = dblInformation + dblIValue;
                if isnan(dblInformation)
                    test = 1;
                    disp('Bad results')
                end
            else
                failcount = 1;
            end
        end
        matAICResults(i,2) = (matAICResults(i,2).*(j-1) - dblInformation + 2*log(dblAIC))./j;
        matAICResults(i,3) = (matAICResults(i,3).*(j-1) + 2*log(dblAIC) - log(10/vecNumOfParts(i)))./j;
        matAICResults(i,4) = (matAICResults(i,4).*(j-1) + 2*log(dblAIC))./j;
    end
    
end
toc
%Plot the figure
plot(log10(matAICResults(:,1)),matAICResults(:,2),'-x',log10(matAICResults(:,1)),matAICResults(:,3),'-o',log10(matAICResults(:,1)),matAICResults(:,4),'-+')
legend('COMIC','AIC $- \ln(|\Omega | / n)$','AIC', 'Interpreter', 'latex')
xlabel('$\log_{10}(n)$','Interpreter', 'latex')
ylabel('Fitness Metric','Interpreter', 'latex')



