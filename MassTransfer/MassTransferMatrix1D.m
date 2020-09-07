function [matMT] = MassTransferMatrix1D(matMT,matT) 
    matMT(:,2) = matT * matMT(:,2);
end

