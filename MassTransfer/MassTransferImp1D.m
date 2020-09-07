function matMT = MassTransferImp1D(matMT, D, dt, celIdx, celRadius)
% first column is position, second is mass
intNumOfPart = size(matMT,1);
%Calculate the local Fickian dispersion (1D), for each particle pair collision
%probability is (rescaled)
% Pij = (1 /(8*pi*D*dt))^(1/2) * exp(-r^2/(8*D*dt))
% Pij = (1/ dblP*pi)^(1/2)*exp(-r^2/dblP), where dblP = 8*D*dt
dblProb = 8*D*dt;
%dblSearchDist = 3*sqrt(dblProb);%this come from Benson code, not sure why
%we choose this to be the search distance. It noted in his code that it is
%for speed.

%Search for nearby paticles with certain radius. We are just the BucketSize
%option for this search, which is a Kd-tree search with maximum of 10 data
%points in the leaf node. There are others search methods.
%[celIdx, celRadius] = rangesearch(matMT(:,1), matMT(:,1), dblSearchDist,'BucketSize',10);

%Loop through all particles and transfer

for i=1:intNumOfPart
    vecNearbyIdx = celIdx{i};
    vecRadius = celRadius{i};
    vecProb = sqrt(1/(dblProb*pi)).*exp(-vecRadius.^2/(dblProb));
    dblProbRescale = sum(vecProb);
    
    %Remove double transfer, by only transfer mass between particles that
    %have position bigger than the current particle.
    vecProb = vecProb(vecNearbyIdx > i );
    vecProb = vecProb./dblProbRescale;
    vecNearbyIdx = vecNearbyIdx(vecNearbyIdx>i);
    intNumOfXPart = length(vecProb);
    if intNumOfXPart > 0
        vecRandIdx = randperm(intNumOfXPart);
        for j=1:intNumOfXPart
            dblMassTrans = 1/2*(matMT(vecNearbyIdx(vecRandIdx(j)),2) - matMT(i,2))*vecProb(vecRandIdx(j));
            %current particle
            matMT(i,2) = max(0,matMT(i,2) + dblMassTrans);
            %the other particle
            matMT(vecNearbyIdx(vecRandIdx(j)),2) = max(0,matMT(vecNearbyIdx(vecRandIdx(j)),2) - dblMassTrans);
        end
        
    end 
end
end


