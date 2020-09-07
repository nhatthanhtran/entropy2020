%Mass transfer D and V recovery: parameters estimation - using the 
%non Gaussian process function
clc; clear all; 
close all;
rng('shuffle')
set(0,'defaultTextInterpreter','latex');
%Initialize
D = 1;
T = 1;
V = 1;
xCenter = 0;
intNumData = 30;%number of data points
intNumOfPart = 3000; %number of particles
intNumOfEns = 1;
Domain = [-6 6];
ZeroDomain = Domain;
dblShift = V*T + xCenter;
Domain = Domain + dblShift;
dx = (Domain(1,2) - Domain(1,1))/(intNumData-1); 

dblBinSize = 0.1;
dblAlpha = 1;
%x = Domain(1,1):dx:Domain(1,2);
%Randomly choosing the x data point within the domain
x = (Domain(1,2) - Domain(1,1))*rand(1,intNumData) + Domain(1,1);
x = 10*rand(1,intNumData) -5;
%x = Domain(1,1):dx:Domain(1,2);
x = sort(x);
intNumData = length(x);

%Get exact solution
vecExactSolution = ExactSolution1D(x,T,D,V,xCenter);


%Multi Gaussian Estimator
intNumOfPart = 3000;
funErrorD = @(d) 1/intNumData*sum(1./(dblAlpha.*vecExactSolution').*(vecExactSolution'-...
    MassTransferImp1DSolution(x,d(1),Domain,ZeroDomain,intNumOfPart,T,intNumOfEns,dblBinSize,d(2),xCenter)).^2);
DValue = [0.5;0.5];
[DValue,DError] = fminsearch(funErrorD,DValue);
DValue
