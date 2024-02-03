% The radar-only system without RIS.
% This is used in the paper: R. Liu, M. Li, Y. Liu, Q. Wu, and Q. Liu, “Joint transmit waveform and passive beamforming design for RIS-aided DFRC systems,”IEEE J. Sel. Topics Signal Process., vol. 16, no .5, pp. 995-1010, Aug. 2022.
% Download this paper at: https://ieeexplore.ieee.org/document/9769997
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-01-28
% Inputs: Prms: the structure of system parameters;
%         Channel: the structure of the channels;
% Outputs: x: transmit waveform; phi: RIS reflection coefficients
%          VSINR: the achieved radar SINR

function [x,VSINR] = get_x_woRIS_radar(Prms,Channel)

M = Prms.M;  L = Prms.L; Q = Prms.Q; sigma2 = Prms.sigma2;
sigmar2 = Prms.sigmar2;  P = Prms.P; Nmax = Prms.Nmax;
res_th = Prms.res_th; ht = Channel.ht;
Hc = Channel.Hc;
Nmax = 100;

A = zeros(M*L,M*L,Q+1);
Jr = zeros(M*L,M*L,Q);
for q = 1:1:Q
    for i = 1:1:M*L
        for j = 1:1:M*L
            if i-j == M*(q-1)
                Jr(i,j,q) = 1;
            end
        end
    end
    A(:,:,q) = kron(eye(L),Hc(q,:)'*Hc(q,:));
end
A(:,:,end) = kron(eye(L),ht'*ht);

%%% initilization
Ht = [Hc;ht];
Ht = Ht/norm(Ht);
xl = get_initial_x_radar(Ht,Prms);
x = sqrt(P/M/L)*kron(ones(L,1),xl);
X = x*x';
Fqphi = zeros(M*L,M*L,Q+1);
Fqphi(:,:,end) = A(:,:,end);
st = Fqphi(:,:,end)*x;
Mt = sigmar2*eye(M*L);
for q = 1:1:Q
    Fqphi(:,:,q) = A(:,:,q)*Jr(:,:,q);
    Mt = Mt + sigma2*Fqphi(:,:,q)*X*Fqphi(:,:,q)';
end
mst = Mt\st;
Mst = mst*mst';

y = x;
mu1 = zeros(M*L,1);
rho = abs(real(st'*mst))/P/10;

Dt = zeros(M*L,M*L);
for q = 1:1:Q
    Dt = Dt + 2*sigma2*Fqphi(:,:,q)'*Mst*Fqphi(:,:,q);
end
Dt = Dt + rho/2*eye(M*L);
dt = -2*Fqphi(:,:,end)'*mst-rho*y + mu1;

Vobj = zeros(1,Nmax);
VSINR = zeros(1,Nmax);
Vres = zeros(1,Nmax);
iter = 1;
res = 1;
while iter <= Nmax && res > res_th

    for i = 1:1:M*L
        Dt(i,i) = real(Dt(i,i)) + 1e-6*real(Dt(1,1));
    end
    sc = norm(dt,2)^2;
    Dt = Dt/sc;
    dt = dt/sc;
    R = chol(Dt);
    %%%% update x
    cvx_begin quiet
    variable x(M*L,1) complex
    minimize real(x'*(R'*R)*x) + real(dt'*x)
    subject to
    abs(x) <= sqrt(P/M/L);
    cvx_end

    %%%% update y
    y = sqrt(P/(M*L))*exp(1i*angle(rho*x+mu1));
    mu1 = mu1 + rho*(x-y);

    X = x*x';
    st = Fqphi(:,:,end)*x;
    Mt = sigmar2*eye(M*L);
    for q = 1:1:Q
        Mt = Mt + sigma2*Fqphi(:,:,q)*X*Fqphi(:,:,q)';
    end
    mst = Mt\st;
    Mst = mst*mst';

    Dt = zeros(M*L,M*L);
    for q = 1:1:Q
        Dt = Dt + 2*sigma2*Fqphi(:,:,q)'*Mst*Fqphi(:,:,q);
    end
    Dt = Dt + rho/2*eye(M*L);
    dt = -2*Fqphi(:,:,end)'*mst-rho*y + mu1;

    Vobj(iter) = real(-st'*mst) + 0.5*rho*(norm(x-y+mu1/rho,2))^2;
    VSINR(iter) = 10*log10(real(sigma2*st'*mst));
    Vres(iter) = norm(x-y,2)^2;
    if iter > 1
        res = abs(1-VSINR(iter)/VSINR(iter-1));
    end
    if iter > 20  && res < 1e-3
        res = abs(1-sum(VSINR(iter-20:1:iter-1))/20/VSINR(iter));
    end
    iter = iter + 1;
end

Vobj(iter:end) = [];
VSINR(iter:end) = [];
Vres(iter:end) = [];





