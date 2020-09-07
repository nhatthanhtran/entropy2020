function [vecApproxSolution] = MassTransferImp1DSolution(vecX,D,vecDomain,vecOrginalDomain,intNumOfPart,T,intNumOfEns,dblBinSize,V,xCenter)
%Function use to create an approximation using the binning method 
%for mass transfer method
% Time to bin
intSizeX = length(vecX);
intNumOfPart = round(intNumOfPart);
vecApproxSolution = zeros(intSizeX,1);
dblCenter = V*T + xCenter;
vecShiftDomain = vecOrginalDomain + dblCenter;
%Binning of thing will affect the velocity approximation
dblBinSize = (vecDomain(1,2) - vecDomain(1,1))/(intNumOfPart-1);

%This is for convenience sake of which method do we want to create our
%numerical concentration. Either work with mass in the mass transfer then
%calculate the concentration later OR directly work with the concentration
strBinningType = 'C'; %M is mass and C is concentration

for j=1:intNumOfEns

    %Use Bin Method to get the approx solution
    matBin = (vecDomain(1,1):dblBinSize:vecDomain(1,2))';
    intBinSize = size(matBin(:,1),1);
    matBin(:,2) = zeros(intBinSize,1);
    vecTempSolutionL = zeros(intSizeX,1);%holding the approx solution for each ensemble
    vecTempSolutionR = zeros(intSizeX,1);
    intNumOfPartAtZero = 1;
    matMTSolution = zeros(intNumOfPart,2);
    %This part need to think about it. This may not the best implementation
    %matMTSolution(:,1) = ((vecDomain(1,2) - vecDomain(1,1))*rand(1,intNumOfPart) + vecDomain(1,1))';
    
    dblSpace = (vecShiftDomain(1,2) - vecShiftDomain(1,1))/(intNumOfPart-1);
    matMTSolution(:,1) = (vecShiftDomain(1,1):dblSpace:vecShiftDomain(1,2))';
    
    %Setting up the initial condition - Direct Delta function 
    %Using first order approx of the initial condition
    if strBinningType == 'M'
        if mod(intNumOfPart,2) == 0
            %Even number of particles, not guarantee the particles at the center
            %Estimate the dirac-delta as first order
            intLeftPart = intNumOfPart/2;
            intRightPart = intLeftPart + 1;
            dblTotalDistance = (dblCenter - matMTSolution(intLeftPart,1))+ ...
                (matMTSolution(intRightPart,1) - dblCenter);
            matMTSolution(intLeftPart,2) = 1 - (dblCenter - matMTSolution(intLeftPart,1))/dblTotalDistance;
            matMTSolution(intRightPart,2) = 1 - (matMTSolution(intRightPart,1) - dblCenter)/dblTotalDistance;
        else%Odd number of particles so we know that their is one at the middle
            matMTSolution(round(intNumOfPart/2),2) = 1;
        end
    elseif strBinningType == 'C'
        
        if mod(intNumOfPart,2) == 0
            intLeftPart = intNumOfPart/2;
            intRightPart = intLeftPart + 1;
            dblTotalDistance = (dblCenter - matMTSolution(intLeftPart,1))+ ...
                (matMTSolution(intRightPart,1) - dblCenter);
            matMTSolution(intLeftPart,2) = (1 - (dblCenter - matMTSolution(intLeftPart,1))/dblTotalDistance)...
                /((matMTSolution(intLeftPart + 1,1) - matMTSolution(intLeftPart - 1,1))/2);
            matMTSolution(intRightPart,2) = (1 - (matMTSolution(intRightPart,1) - dblCenter)/dblTotalDistance)...
                /((matMTSolution(intRightPart + 1,1) - matMTSolution(intRightPart - 1,1))/2);
        else
            intCenterPart = round(intNumOfPart/2);
            %calculate the concetration of particle at the center
            matMTSolution(intCenterPart,2) = 1/ ((matMTSolution(intCenterPart + 1,1) - matMTSolution(intCenterPart - 1,1))/2);
            
        end
    else
        disp('Bad bin type input')
        break;
    end
    
    %Start the mass/concentration transfer
    intNumOfTimeSteps = 100;
    dt = T/(intNumOfTimeSteps);
    if ((vecDomain(1,2) - vecDomain(1,1))/intNumOfPart) > sqrt(2*D*dt)
        disp('Warning: Bad particles spacing.')
    end
    % Figure out the nearby particles
    dblProb = 8*D*dt;
    dblSearchDist = 3*sqrt(dblProb);
    %dblSearchDist = 20;
    %Search for nearby paticles with certain radius. We are just the BucketSize
    %option for this search, which is a Kd-tree search with maximum of 10 data
    %points in the leaf node. There are others search methods.
    [celIdx, celRadius] = rangesearch(matMTSolution(:,1), matMTSolution(:,1), dblSearchDist,'BucketSize',max(10,round(intNumOfPart*1e-2)));
    %[celIdx, celRadius] = rangesearch(matMTSolution(:,1), matMTSolution(:,1), dblSearchDist,'BucketSize',100);
    matT = MassTransferProbMat1D(matMTSolution,D,dt,celIdx,celRadius);
    for k = 1:intNumOfTimeSteps
        %matMTSolution = MassTransferImp1D(matMTSolution,D,dt,celIdx,celRadius);
        %matMTSolution = MassTransferMatrix1D(matMTSolution,D,dt,celIdx,celRadius);
        matMTSolution = MassTransferMatrix1D(matMTSolution,matT);
    end
%     vecApproxSolution = matMTSolution(:,2);
%     break;
    %Compute the approximation
    if strBinningType == 'C'
        vecApproxTemp = MassTransferBinning1D(vecX,matMTSolution);
        vecApproxSolution = (vecApproxSolution.*(j-1) + vecApproxTemp)./j;
    elseif strBinningType == 'M'
        %     %Left bin [ )
        intXIdx = 1;
        for i=1:intBinSize-1
            vecAboveMass = matMTSolution(matMTSolution(:,1) >= matBin(i,1),2);
            vecBelowMassP1 = matMTSolution(matMTSolution(:,1) < matBin(i+1,1),2);
            matBin(i,2) = (sum(vecAboveMass) + ...
                sum(vecBelowMassP1)-1)/dblBinSize;
            
            %Left bin approx - getting all of the approx solution that land in certain bin
            while vecX(intXIdx) >= matBin(i,1) && vecX(intXIdx) < matBin(i+1,1)
                vecTempSolutionL(intXIdx,1) = matBin(i,2);
                intXIdx = intXIdx + 1;
                if intXIdx > intSizeX %there no need to keep binning as we get the approx for all x points
                    break
                end
            end
            if intXIdx > intSizeX
                break
            end
        end
        
        
        %Catch the last xindex if it on the boundary
        if intXIdx <= intSizeX
            vecAboveMass = matMTSolution(matMTSolution(:,1) >= matBin(intBinSize,1),2);
            vecBelowMassP1 = matMTSolution(matMTSolution(:,1) < (matBin(intBinSize,1) + dblBinSize),2);
            dblLastBinValue = (sum(vecAboveMass) + ...
                sum(vecBelowMassP1)-1)/dblBinSize;
            while intXIdx <= intSizeX
                vecTempSolutionL(intXIdx,1) = dblLastBinValue;
                intXIdx = intXIdx + 1;
            end
        end
        
        %Do the right bin ( ]
        intXIdx = 1;
        %Get the first right bound
        vecAboveMass = matMTSolution(matMTSolution(:,1) > (matBin(1,1) - dblBinSize),2);
        vecBelowMassP1 = matMTSolution(matMTSolution(:,1) <= matBin(1,1),2);
        dblFirstBinValue = (sum(vecAboveMass) + ...
            sum(vecBelowMassP1)-1)/dblBinSize;
        while vecX(intXIdx) == vecDomain(1,1)
            vecTempSolutionR(intXIdx,1) = dblFirstBinValue;
            intXIdx = intXIdx + 1;
        end
        
        %Get the rest of the bin
        for i=1:intBinSize-1
            vecAboveMass = matMTSolution(matMTSolution(:,1) > matBin(i,1),2);
            vecBelowMassP1 = matMTSolution(matMTSolution(:,1) <= matBin(i+1,1),2);
            matBin(i,2) = (sum(vecAboveMass) + ...
                sum(vecBelowMassP1)-1)/dblBinSize;
            
            %Left bin approx - getting all of the approx solution that land in certain bin
            while vecX(intXIdx) > matBin(i,1) && vecX(intXIdx) <= matBin(i+1,1)
                vecTempSolutionR(intXIdx,1) = matBin(i,2);
                intXIdx = intXIdx + 1;
                if intXIdx > intSizeX %there no need to keep binning as we get the approx for all x points
                    break
                end
            end
            if intXIdx > intSizeX
                break
            end
        end
        
        vecApproxSolution = (vecApproxSolution.*(j-1) + (vecTempSolutionL + vecTempSolutionR)./2)./j;
    else
        disp('Bad bin type input')
        break;
    end
end
end

