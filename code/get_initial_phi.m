% The initialization for phi by solving (62).
% This is used in the paper: R. Liu, M. Li, Y. Liu, Q. Wu, and Q. Liu, “Joint transmit waveform and passive beamforming design for RIS-aided DFRC systems,”IEEE J. Sel. Topics Signal Process., vol. 16, no .5, pp. 995-1010, Aug. 2022.
% Download this paper at: https://ieeexplore.ieee.org/document/9769997
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-01-28
% Inputs: Prms: the structure of system parameters;%
% Outputs: phi: RIS reflection coefficients
function phi = get_initial_phi(Prms,B1,C1,B2,C2)
N = Prms.N;  Q = Prms.Q;

% Create the problem structure.
manifold = complexcirclefactory(N);
problem.M = manifold;
% Define the problem cost function and its gradient.

problem.cost = @cost;
    function f = cost(v)
        f = v'*C1*v + 2*real(B1*v);
        for q = 1:1:Q
            f = f - v'*C2(:,:,q)*v - 2*real(B2(q,:)*v);
        end

    end
problem.grad = @(v) problem.M.egrad2rgrad(v,egrad(v));
    function g = egrad(v)
        g = 2*C1*v + 2*B1';
        for q = 1:1:Q
            g = g - 2*C2(:,:,q)*v - 2*B2(q,:)';
        end
    end

% Execute the optimization
options.tolgradnorm = 1e-3;
options.maxiter = 1000;
options.minstepsize = 1e-3;
options.verbosity = 0;

[v,aa,bb] = conjugategradient(problem,[],options);

phi = conj(v);
end


