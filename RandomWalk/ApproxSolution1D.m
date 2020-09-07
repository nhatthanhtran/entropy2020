function [vecApproxSolution] = ApproxSolution1D(vecX,dblBinSize,vecDomain,D,T,intNumOfPart,intNumEns,V,xCenter,strBinType)
%return the approximation of the concentration
%This is for random walk method
%BinType options: 'L' - [) or 'R' - (]
intSizeX = length(vecX);
 
intNumOfPart = round(intNumOfPart);% this is because minsearch may not give 
%an integer value
vecApproxSolution = zeros(intSizeX,1);
vecApproxSolutionL = zeros(intSizeX,1);
vecApproxSolutionR = zeros(intSizeX,1);
if strBinType == 'L'|| strBinType == 'B'%left bin type
for j=1:intNumEns

    %Use Bin Method to get the approx solution
    matBin = (vecDomain(1,1):dblBinSize:vecDomain(1,2))';
    intBinSize = size(matBin(:,1),1);
    matBin(:,2) = zeros(intBinSize,1);
    vecTempSolution = zeros(intSizeX,1);%holding the approx solution for each ensemble
    vecRWSolution = RandomWalk1D(intNumOfPart,D,T,V,xCenter);
    intXIdx = 1;
    for i=1:intBinSize-1
       matBin(i,2) = ((sum(vecRWSolution>=matBin(i,1)) + sum(vecRWSolution<matBin(i+1,1))-intNumOfPart)...
           /intNumOfPart)/dblBinSize;

       %Left bin approx - getting all of the approx solution that land in certain bin 
       while vecX(intXIdx) >= matBin(i,1) && vecX(intXIdx) < matBin(i+1,1)
           vecTempSolution(intXIdx,1) = matBin(i,2);
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
        dblLastBinValue = ((sum(vecRWSolution>=matBin(intBinSize,1)) + ...
            sum(vecRWSolution<(matBin(intBinSize,1)+dblBinSize)) - intNumOfPart)...
           /intNumOfPart)/dblBinSize;
       while intXIdx <= intSizeX
          vecTempSolution(intXIdx,1) = dblLastBinValue;
          intXIdx = intXIdx + 1;
       end
    end
    
    %matBin(intBinSize,2) = matBin(intBinSize-1,2);

    vecApproxSolutionL = (vecApproxSolutionL.*(j-1) + vecTempSolution)./j;
    %Other methods beside binning can be use - later
end
end
%try right bin method
if strBinType == 'R' || strBinType == 'B'
for j=1:intNumEns

    %Use Bin Method to get the approx solution
    matBin = (vecDomain(1,1):dblBinSize:vecDomain(1,2))';
    intBinSize = size(matBin(:,1),1);
    matBin(:,2) = zeros(intBinSize,1);
    vecTempSolution = zeros(intSizeX,1);%holding the approx solution for each ensemble
    vecRWSolution = RandomWalk1D(intNumOfPart,D,T,V,xCenter);
    intXIdx = 1;
    
    %Get the first left bound
    dblFirstBinValue = ((sum(vecRWSolution>(matBin(1,1)-dblBinSize)) + ...
    sum(vecRWSolution<=(matBin(1,1))) - intNumOfPart)...
    /intNumOfPart)/dblBinSize;
    while vecX(intXIdx) == vecDomain(1,1)
        vecTempSolution(intXIdx,1) = dblFirstBinValue;
        intXIdx = intXIdx + 1;
    end

    for i=1:intBinSize-1
       matBin(i,2) = ((sum(vecRWSolution>matBin(i,1)) + sum(vecRWSolution<=matBin(i+1,1))-intNumOfPart)...
           /intNumOfPart)/dblBinSize;

       %Left bin approx - getting all of the approx solution that land in certain bin 
       while vecX(intXIdx) > matBin(i,1) && vecX(intXIdx) <= matBin(i+1,1)
           vecTempSolution(intXIdx,1) = matBin(i,2);
           intXIdx = intXIdx + 1;
           if intXIdx > intSizeX %there no need to keep binning as we get the approx for all x points
               break
           end
       end
       if intXIdx > intSizeX
          break 
       end
    end

    vecApproxSolutionR = (vecApproxSolutionR.*(j-1) + vecTempSolution)./j;
    %Other methods beside binning can be use - later
end
end

%Get the result back
if strBinType == 'L'
    vecApproxSolution = vecApproxSolutionL; 
elseif strBinType == 'R'
    vecApproxSolution = vecApproxSolutionR;
elseif strBinType == 'B'
    vecApproxSolution = (vecApproxSolutionL + vecApproxSolutionR)./2;
else
    disp('Incorrect BinType provided.')
end
end

