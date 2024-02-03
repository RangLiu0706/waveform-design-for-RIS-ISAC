% The initialization for x.
% This is used in the paper: R. Liu, M. Li, Y. Liu, Q. Wu, and Q. Liu, “Joint transmit waveform and passive beamforming design for RIS-aided DFRC systems,”IEEE J. Sel. Topics Signal Process., vol. 16, no .5, pp. 995-1010, Aug. 2022.
% Download this paper at: https://ieeexplore.ieee.org/document/9769997
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-01-28
% Inputs: Prms: the structure of system parameters;%
% Outputs: x: radar waveform
function x = get_initial_x_radar(Ht,Prms)
M = Prms.M;  Q = Prms.Q;

% Create the problem structure.
manifold = complexcirclefactory(M);
problem.M = manifold;
% Define the problem cost function and its gradient.

problem.cost = @cost;
    function f = cost(x)
        f = x'*Ht(end,:)'*Ht(end,:)*x;
        for q = 1:1:Q
            f = f - x'*Ht(q,:)'*Ht(q,:)*x;
        end

    end
problem.grad = @(x) problem.M.egrad2rgrad(x,egrad(x));
    function g = egrad(x)
        g = 2*Ht(end,:)'*Ht(end,:)*x;
        for q = 1:1:Q
            g = g - 2*Ht(q,:)'*Ht(q,:)*x;
        end
    end

% Execute the optimization
options.tolgradnorm = 1e-6;
options.maxiter = 1000;
options.minstepsize = 1e-6;
options.verbosity = 0;

[x,~,~] = conjugategradient(problem,[],options);

end


