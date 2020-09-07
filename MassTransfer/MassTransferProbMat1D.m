function [matT] = MassTransferProbMat1D(matMT, D, dt, celIdx, celRadius)
%Using Mike code with some modification
    N = length(matMT(:,1));
    beta = 1;
    dblProb = (beta^(-1))*4*D*dt;%This is an update delete the line 5 if want to go back and beta = 1/2
    Nclose = sum(cellfun('length', celIdx));
    row = zeros(1, Nclose);
    col = zeros(1, Nclose);
    val = zeros(1, Nclose);

    start = 1;

    for ii = 1 : N
        finish = start - 1 + length(celIdx{ii});
        row(start : finish) = ii;
        col(start : finish) = celIdx{ii};
        val(start : finish) = celRadius{ii};
        start = finish + 1;
    end
    
    clear celIdx celRadius
     
%     exponentiate before sparsing to make life easier
%sqrt(1/(dblProb*pi)).*exp(-vecRadius.^2/(dblProb));
    %val =  factor .* exp(-(val.^2) ./ denom);  
    val = sqrt(1/(dblProb*pi)).*exp(-(val.^2)/(dblProb));
    Pmat = sparse(row, col, val);
%    clear row col val
    clear val
    
%     normalize Pmat such that it is symmetric, using the arithmetic mean of
%     colsum and rowsum
    colsum = sparse(sum(Pmat));
    rowsum = sparse(sum(Pmat, 2));
    
% %     this is the simpler, memory-intensive way of doing things, and
% %     builds a fully-dense normMat
%     normMat = (colsum + rowsum) ./ 2;
%     
%     clear colsum rowsum
%     
%     normMat(normMat < 1e-16) = 0;
%     normMat(normMat == 0) = 1;
%     Pmat = Pmat ./ normMat;
    
%     %     this requires some extra ./'s to avoid divsion by zero and
%     %     de-sparsing things, but it keeps it memory efficient
%     Pmat = spfun(@(x) 1 ./ x, (2 .* Pmat) ./ rowsum) + spfun(@(x) 1 ./ x, (2 .* Pmat) ./ colsum);
% %     Pmat = 1 ./ ((2 .* Pmat) ./ rowsum) + 1 ./ ((2 .* Pmat) ./ colsum);
%     
% %     
%     clear rowsum colsum

    if N > 1e5
        rowNorm = spfun(@(x) 1 ./ x, (2 .* Pmat) ./ rowsum);
        colNorm = spfun(@(x) 1 ./ x, (2 .* Pmat) ./ colsum);
        for i = 1 : Nclose
            Pmat(row(i), col(i)) = rowNorm(row(i), col(i)) + colNorm(row(i), col(i));
        end
        clear rowNorm colNorm
    else
        Pmat = spfun(@(x) 1 ./ x, (2 .* Pmat) ./ rowsum) + spfun(@(x) 1 ./ x, (2 .* Pmat) ./ colsum);
    end

    clear row col

    clear rowsum colsum

    
    Pmat = spfun(@(x) 1 ./ x, Pmat);
    
%     This is the generalization of the explicit mass-transfer matrix from
%     the SPH paper from Guillem and Mike
    Tmat = speye(N) - beta * spdiags(Pmat * ones(N, 1), 0, N, N) + beta * Pmat;
    matT = Tmat;
    clear Tmat;
    clear Pmat
end

