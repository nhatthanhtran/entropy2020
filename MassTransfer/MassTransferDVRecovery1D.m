%Mass transfer D and V recovery: parameters estimation
clc; clear all; 
close all;
rng('shuffle')
set(0,'defaultTextInterpreter','latex');
%Initialize
D = 1;
T = 1;
V = 1;
xCenter = 0;
intNumData = 6;%number of data points
intNumOfPart = 3000; %number of particles
intNumOfEns = 1;
Domain = [-6 6];
ZeroDomain = Domain;
dblShift = V*T + xCenter;
Domain = Domain + dblShift;
dx = (Domain(1,2) - Domain(1,1))/(intNumData-1); 

bin = 0.1;%dummy variable to run the code, wont do anything - legacy paramter
%x = Domain(1,1):dx:Domain(1,2);
%Randomly choosing the x data point within the domain
%Can choose between equally space or random
%x = (Domain(1,2) - Domain(1,1))*rand(1,intNumData) + Domain(1,1);
%x = 10*rand(1,intNumData) - 5;
%x = Domain(1,1):dx:Domain(1,2);
x = -5:10/(intNumData-1):5;
x = sort(x);

intNumData = length(x);
%Get exact solution
vecExactSolution = ExactSolution1D(x,T,D,V,xCenter);

%Find D and V
funErrorDV = @(dv) 1/intNumData*norm(vecExactSolution' -...
    MassTransferImp1DSolution(x,dv(1),Domain,ZeroDomain,intNumOfPart,T,intNumOfEns,bin,dv(2),xCenter),2)^2;
DVValue = [0.5;0.5];
[DVValue,DVError] = fminsearch(funErrorDV,DVValue);
format long e
DVValue


