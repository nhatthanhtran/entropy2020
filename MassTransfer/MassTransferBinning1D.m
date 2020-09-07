function [vecApproxSolution] = MassTransferBinning1D(vecX,matConcentration)
%Use to work with concetration and compute the concentration at certain
%spatial position
%Bin the solution
vecApproxSolution = zeros(length(vecX),1);
intOrderOfApprox = 1;%Control the order of approximation that use to
%compute the concentration

%TRY zero order
if intOrderOfApprox == 0
    vecSupportBin = zeros(length(matConcentration(:,1))+1,1);
    vecSupportBin(1,1) = matConcentration(1,1);
    vecSupportBin(end,1) = matConcentration(end,1);
    for q = 1:length(matConcentration(:,1))-1
        vecSupportBin(q+1,1) = (matConcentration(q+1,1) + matConcentration(q,1))/2;
    end
    intXIdx = 1;
    intBinIdx = 1;
    while intXIdx <= length(vecX)
        blnBinFound = 0;
        while ~blnBinFound
            if vecX(intXIdx) >= vecSupportBin(intBinIdx) && vecX(intXIdx) <= vecSupportBin(intBinIdx+1)
                vecApproxSolution(intXIdx) = matConcentration(intBinIdx,2);
                blnBinFound = 1;
                intXIdx = intXIdx + 1;
            else
                intBinIdx = intBinIdx + 1;
            end
        end
    end
    
elseif intOrderOfApprox == 1
    % TRY first order approx
    % %Particles outside of the domain assign to be zero?
    vecLeftIdx = vecX < matConcentration(1,1);
    vecRightIdx = vecX > matConcentration(end,1);
    intXIdx = 0;
    if sum(vecLeftIdx) == 0
        intXIdx = 1;
    else
        intXIdx = find(vecLeftIdx, 1, 'last')+1;
    end
    
    intEndIdx = 0;
    if sum(vecRightIdx) == 0
        intEndIdx = length(vecRightIdx);
    else
        intEndIdx = find(vecRightIdx, 1, 'first')-1;
    end
    
    intBinIdx = 1;
    while intXIdx <= intEndIdx
        blnBinFound = 0;
        while ~blnBinFound
            %the position is exactly the particle
            if vecX(intXIdx) == matConcentration(intBinIdx,1)
                vecApproxSolution(intXIdx) = matConcentration(intBinIdx,2);
                intXIdx = intXIdx + 1;
                blnBinFound = 1;
                %the position is in between the particles
            elseif vecX(intXIdx) >= matConcentration(intBinIdx,1) && vecX(intXIdx) <= matConcentration(intBinIdx+1,1)
                dblDistance = matConcentration(intBinIdx+1,1) - matConcentration(intBinIdx,1);
                %first order approx - weighted average base on distance
                vecApproxSolution(intXIdx) = (1-(vecX(intXIdx) - matConcentration(intBinIdx,1))./dblDistance).*matConcentration(intBinIdx,2) +...
                    (1-(matConcentration(intBinIdx+1,1) - vecX(intXIdx))./dblDistance).*matConcentration(intBinIdx + 1,2);
                intXIdx = intXIdx + 1;
                blnBinFound  = 1;
            else
                intBinIdx = intBinIdx + 1;
            end
            
            if intBinIdx > length(matConcentration(:,1))
                disp("Can't find all the bin for the x positions - Break Infinite Loop");
                disp("Fail at index - " + num2str(intXIdx));
                break;
            end
        end
        if intBinIdx > length(matConcentration(:,1))
            break;
        end
    end
else
    disp('Need to config the order of approximation')
end
end

