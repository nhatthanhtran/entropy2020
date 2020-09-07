%Parameter estimation for random walk
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
intNumData = 30;%number of data points - need to modify this for different number of data
intNumOfPartStar = 5000;%NEED to change this according to the data value
intNumOfEns = 15; %number of ensemble for random walk
Domain = [-5 5] + dblXShift;
%dx = (Domain(1,2) - Domain(1,1))/(intNumData-1); 
dblBinSize = 0.1;

%Randomly choosing the x data point within the domain
x = (Domain(1,2) - Domain(1,1))*rand(1,intNumData) + Domain(1,1);
%x = Domain(1,1):dx:Domain(1,2);
x = sort(x);

intNumData = length(x);

%Get exact solution
vecExactSolution = ExactSolution1D(x,T,D,V,xCenter);

intNumOfTrials = 30;

matResults = zeros(intNumOfTrials,4);

for i=1:intNumOfTrials
    
    funErrorDStar =@(dv) 1/intNumData*norm(vecExactSolution'-...
        ApproxSolution1D(x,dblBinSize,Domain,dv(1),T,intNumOfPartStar,intNumOfEns,dv(2),xCenter,'B'),2)^2;
    DValueStar = [0.5;0.5];
    
    [DValueStar, DErrorStar] = fminsearch(funErrorDStar,DValueStar);
    matResults(i,1) = DValueStar(1);
    matResults(i,3) = DValueStar(2);
end

figure(1)
subplot(1,2,1)
boxplot(matResults(:,1),'Labels',{''})

ylabel('Estimate $D$','Interpreter', 'latex')
subplot(1,2,2)
boxplot(matResults(:,3),'Labels',{''})
ylabel('Estimate $v$','Interpreter', 'latex')

