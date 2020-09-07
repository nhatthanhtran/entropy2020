function [vecPosition] = RandomWalk1D(N, D, T, V, xCenter)
%return final position of particles
%Time step control, help with control the error of the randomness
intNumOfSteps = 25;
dt = T/intNumOfSteps;
%dt = 0.02;
vecPosition = zeros(N,1) + xCenter;

for i=1:intNumOfSteps
    vecMove = sqrt(2*D*dt).*randn(N,1);
    vecPosition = vecPosition + vecMove;
end
vecPosition = vecPosition + V*T;
end

