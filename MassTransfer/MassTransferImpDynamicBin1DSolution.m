function [vecApproxSolution,matBin,vecMass,intCenterIdx] = MassTransferImpDynamicBin1DSolution(vecX,D,vecDomain,vecOrginalDomain,intNumOfPart,T,intNumOfEns,vecBin,V,xCenter)
% Use to create the approx solution for mass transfer with dynamic binning
intSizeX = length(vecX);
intNumOfPart = round(intNumOfPart);
vecApproxSolution = zeros(intSizeX,1);
dblCenter = V*T + xCenter;
vecShiftDomain = vecOrginalDomain + dblCenter;
%Binning of thing will affect the velocity approximation
dblBinSize = (vecDomain(1,2) - vecDomain(1,1))/(intNumOfPart);
for j=1:intNumOfEns
    
    %Use Bin Method to get the approx solution
    matBin = vecBin';
    intBinSize = size(matBin(:,1),1);
    matBin(:,2) = zeros(intBinSize,1);
    intNumOfPartAtZero = 1;
    matMTSolution = zeros(intNumOfPart,2);
    dblSpace = (vecShiftDomain(1,2) - vecShiftDomain(1,1))/(intNumOfPart-1);
    
    matMTSolution(:,1) = vecBin';
    
    %Set up the concentration
    for z = 1:length(matMTSolution(:,1))
        if matMTSolution(z,1) > dblCenter
            intLeftPart = z-1;
            intRightPart = z;
            dblTotalDistance = (dblCenter - matMTSolution(intLeftPart,1))+ ...
                (matMTSolution(intRightPart,1) - dblCenter);
            
            dblLSupportVolume = (matMTSolution(z-1,1) - matMTSolution(z-2,1))/2 + (matMTSolution(z,1) - matMTSolution(z-1,1))/2;
            dblRSupportVolume = (matMTSolution(z,1) - matMTSolution(z-1,1))/2 + (matMTSolution(z+1,1) - matMTSolution(z,1))/2;
            matMTSolution(intLeftPart,2) = (1 - (dblCenter - matMTSolution(intLeftPart,1))/dblTotalDistance)/dblLSupportVolume;
            matMTSolution(intRightPart,2) = (1 - (matMTSolution(intRightPart,1) - dblCenter)/dblTotalDistance)/dblRSupportVolume;
            intCenterIdx = z;
            break;
        end
    end
    
    intNumOfTimeSteps = 100;
    dt = T/intNumOfTimeSteps;
    if ((vecDomain(1,2) - vecDomain(1,1))/intNumOfPart) > sqrt(2*D*dt)
        disp('Warning: Bad particles spacing.')
    end
    % Figure out the nearby particles
    dblProb = 8*D*dt;
    dblSearchDist = 3*sqrt(dblProb);
    %dblSearchDist = 1;
    %Search for nearby paticles with certain radius. We are just the BucketSize
    %option for this search, which is a Kd-tree search with maximum of 10 data
    %points in the leaf node. There are others search methods.
    [celIdx, celRadius] = rangesearch(matMTSolution(:,1), matMTSolution(:,1), dblSearchDist,'BucketSize',max(10,round(intNumOfPart*1e-2)));
    matT = MassTransferProbMat1D(matMTSolution,D,dt,celIdx,celRadius);
    for k = 1:intNumOfTimeSteps
        %matMTSolution = MassTransferImp1D(matMTSolution,D,dt,celIdx,celRadius);
        %matMTSolution = MassTransferMatrix1D(matMTSolution,D,dt,celIdx,celRadius);
        matMTSolution = MassTransferMatrix1D(matMTSolution,matT);
    end
    vecMass = matMTSolution(:,2);
    
    %Dont need left or right bin?
    
    %Create a support volume bin
    vecTempApproxSolution = zeros(intSizeX,1);
    vecTempApproxSolution = MassTransferBinning1D(vecX,matMTSolution);
%     vecSupportBin = zeros(intNumOfPart+1,1);
%     vecSupportBin(1,1) = matMTSolution(1,1);
%     vecSupportBin(end,1) = matMTSolution(end,1);
%     for q = 1:intNumOfPart-1
%         vecSupportBin(q+1,1) = (matMTSolution(q+1,1) + matMTSolution(q,1))/2;
%     end
%     intXIdx = 1;
%     intBinIdx = 1;
%     while intXIdx <= length(vecX)
%         blnBinFound = 0;
%         while ~blnBinFound
%             if vecX(intXIdx) >= vecSupportBin(intBinIdx) && vecX(intXIdx) <= vecSupportBin(intBinIdx+1)
%                 vecTempApproxSolution(intXIdx) = matMTSolution(intBinIdx,2);
%                 blnBinFound = 1;
%                 intXIdx = intXIdx + 1;
%             else
%                 intBinIdx = intBinIdx + 1;
%             end
%         end
%         
%         
%     end
    vecApproxSolution = (vecApproxSolution.*(j-1) + vecTempApproxSolution)./j;
end
end

