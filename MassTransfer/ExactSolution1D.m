function vecExactSolution = ExactSolution1D(vecX, T, D, V, xCenter)
%exact solution of the 1D diffusion eq with dirac delta initial condition
    funExactSolution = @(x) 1/sqrt(4*pi*D*T).*exp(-(abs(x - V*T - xCenter).^2)/(4*D*T));
    vecExactSolution = feval(funExactSolution, vecX);
end

