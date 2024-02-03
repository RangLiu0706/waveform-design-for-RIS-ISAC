% The proposed joint transmit waveform and passive beamforming design algorithm.
% This is used in the paper: R. Liu, M. Li, Y. Liu, Q. Wu, and Q. Liu, “Joint transmit waveform and passive beamforming design for RIS-aided DFRC systems,”IEEE J. Sel. Topics Signal Process., vol. 16, no .5, pp. 995-1010, Aug. 2022.
% Download this paper at: https://ieeexplore.ieee.org/document/9769997
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-01-28
% Inputs: Prms: the structure of system parameters;
%         Channel: the structure of the channels;
%         S: the transmitted communication symbols
% Outputs: x: transmit waveform; phi: RIS reflection coefficients
%          VSINR: the achieved radar SINR

function [x,phi,VSINR] = get_x_phi_CI(Prms,Channel,S)

M = Prms.M; N = Prms.N; K = Prms.K; L = Prms.L; Q = Prms.Q; sigma2 = Prms.sigma2;
sigmar2 = Prms.sigmar2; Phi = Prms.Phi; P = Prms.P; Nmax = Prms.Nmax;
res_th = Prms.res_th; gamma = Prms.gamma; ht = Channel.ht; hrt = Channel.hrt;
Hc = Channel.Hc; Hrc = Channel.Hrc; G = Channel.G; Hu = Channel.Hu;
Hru = Channel.Hru;

U = [eye(N)  1i*eye(N)];
A = zeros(M*L,M*L,Q+1);
B = zeros(N,M,Q+1);
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
    B(:,:,q) = diag(Hrc(q,:))*G;
end
A(:,:,end) = kron(eye(L),ht'*ht);
B(:,:,end) = diag(hrt)*G;
Iell = eye(L);
Gk = zeros(N*L,M*L,K);
for k = 1:1:K
    Gk(:,:,k) = kron(eye(L),diag(Hru(k,:))*G);
end
Hkl = zeros(K*L,M*L);
gamma_CI = zeros(K*L,2);
for ell = 1:1:L
    for k = 1:1:K
        Hkl((ell-1)*K+k,:) = kron(Iell(ell,:),Hu(k,:));
        gamma_CI((ell-1)*K+k,1) = exp(-1i*angle(S(k,ell)))*(sin(Phi)+exp(-1i*pi/2)*cos(Phi))/(gamma(k)*sin(Phi));
        gamma_CI((ell-1)*K+k,2) = exp(-1i*angle(S(k,ell)))*(sin(Phi)-exp(-1i*pi/2)*cos(Phi))/(gamma(k)*sin(Phi));
    end
end
sc = abs(Hkl(1,1));
Hkl = Hkl/sc;
Gk = Gk/sc;
gamma_CI = gamma_CI*sc;

%%% initilize phi
B2 = zeros(Q,N);
C2 = zeros(N,N,Q);
B1 = ht*G'*diag(hrt');
C1 = diag(hrt)*(G*G')*diag(hrt');
for q = 1:1:Q
    B2(q,:) = Hc(q,:)*G'*diag(Hrc(q,:)');
    C2(:,:,q) = diag(Hrc(q,:))*(G*G')*diag(Hrc(q,:)');
end
sc = norm(B1);
B1 = B1./sc;
C1 = C1./sc;
B2 = B2./sc;
C2 = C2./sc;
phi = get_initial_phi(Prms,B1,C1,B2,C2);
Phi2 = phi*phi.';
phi_bar = [real(phi);imag(phi)];
Hkl_tilde = zeros(K*L,M*L);
for ell = 1:1:L
    for k = 1:1:K
        Hkl_tilde((ell-1)*K+k,:) = Hkl((ell-1)*K+k,:) + kron(Iell(ell,:),phi.')*Gk(:,:,k);
    end
end
%%% initial x
Ht = [diag(gamma_CI(:,1))*Hkl_tilde;diag(gamma_CI(:,2))*Hkl_tilde];
Ht = Ht./max(max(abs(Ht)));

% x = get_initial_x(Ht);
cvx_begin quiet
variable x(M*L,1) complex
variable t
maximize t
subject to
t <= real(Ht*x);
abs(x) <= sqrt(P/M/L);
cvx_end
% x = sqrt(P/M/L)*x;
X = x*x';
Fqphi = zeros(M*L,M*L,Q+1);
Fqphi(:,:,end) = A(:,:,end) + kron(eye(L),ht'*phi.'*B(:,:,end)) + kron(eye(L),B(:,:,end)'*phi*ht)...
    + kron(eye(L),B(:,:,end)'*Phi2*B(:,:,end));
st = Fqphi(:,:,end)*x;
Mt = sigmar2*eye(M*L);
for q = 1:1:Q
    Fqphi(:,:,q) = ( A(:,:,q) + kron(eye(L),Hc(q,:)'*phi.'*B(:,:,q)) + kron(eye(L),B(:,:,q)'*phi*Hc(q,:))...
        + kron(eye(L),B(:,:,q)'*Phi2*B(:,:,q)) )*Jr(:,:,q);
    Mt = Mt + sigma2*Fqphi(:,:,q)*X*Fqphi(:,:,q)';
end
mst = Mt\st;
Mst = mst*mst';

z = phi;
y = x;
mu1 = zeros(M*L,1);
mu2 = zeros(N,1);
rho = abs(real(st'*mst))/(P+N*L);

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
    real(gamma_CI(:,1).*(Hkl_tilde*x)) >= 1;
    real(gamma_CI(:,2).*(Hkl_tilde*x)) >= 1;
    abs(x) <= sqrt(P/M/L);
    cvx_end

    X = x*x';
    st = Fqphi(:,:,end)*x;
    Mt = sigmar2*eye(M*L);
    for q = 1:1:Q
        Fqphi(:,:,q) = ( A(:,:,q) + kron(eye(L),Hc(q,:)'*phi.'*B(:,:,q)) + kron(eye(L),B(:,:,q)'*phi*Hc(q,:))...
            + kron(eye(L),B(:,:,q)'*Phi2*B(:,:,q)) )*Jr(:,:,q);
        Mt = Mt + sigma2*Fqphi(:,:,q)*X*Fqphi(:,:,q)';
    end
    mst = Mt\st;
    Gkl_tilde = zeros(K*L,N);
    for ell = 1:1:L
        for k = 1:1:K
            Gkl_tilde((ell-1)*K+k,:) = (reshape(Gk(:,:,k)*x,N,L)*Iell(:,ell)).';
        end
    end
    %%%%% update phi
    X0 = reshape(x,M,L);

    Rmst = reshape(mst,M,L);
    ft = -2*(conj(B(:,:,end)*X0*Rmst'*ht') + B(:,:,end)*Rmst*X0'*ht');
    Ft = zeros(N,N);
    Ftv_tilde = -2*B(:,:,end)*Rmst*X0'*B(:,:,end)';
    lambda1 = 0;
    for q = 1:1:Q
        aq = A(:,:,q)*Jr(:,:,q)*x;
        Xq = reshape(Jr(:,:,q)*x,M,L);
        temp1 = conj(B(:,:,q)*Xq*Rmst'*Hc(q,:)') + B(:,:,q)*Rmst*Xq'*Hc(q,:)';
        Temp2 = B(:,:,q)*Rmst*Xq'*B(:,:,q)';
        ft = ft + 4*sigma2*temp1*mst'*aq;
        Ft = Ft + 2*sigma2*(temp1*temp1');
        Ftv_tilde = Ftv_tilde + 4*sigma2*mst'*aq*Temp2 + 4*sigma2*phi.'*Temp2'*phi*Temp2;
        lambda1 = lambda1 + 2*sigma2*norm(Temp2,'fro')^2;
    end
    Ftv_tilde = Ftv_tilde -2*lambda1*Phi2;
    Ftv_bar = [real(Ftv_tilde) imag(Ftv_tilde);imag(Ftv_tilde) -real(Ftv_tilde)];
    [~,eD1] = eigsort(Ftv_bar+Ftv_bar.');
    lambda2 = max(diag(eD1));
    ftv_bar = (Ftv_bar + Ftv_bar.' - lambda2*eye(2*N))*phi_bar;

    gradf = zeros(2*N,1);
    Hessf = zeros(2*N,2*N);
    for q = 1:1:Q
        Xq = reshape(Jr(:,:,q)*x,M,L);
        cq = 4*sigma2*(conj(B(:,:,q)*Xq*Rmst'*Hc(q,:)') + B(:,:,q)*Rmst*Xq'*Hc(q,:)');
        Bq = B(:,:,q)*Xq*Rmst'*B(:,:,q)';
        cq1 = [real(cq);imag(cq)];
        cq2 = [imag(cq);-real(cq)];
        Bq1 = [real(Bq) imag(Bq);imag(Bq) -real(Bq)];
        Bq2 = [imag(Bq) -real(Bq);-real(Bq) -imag(Bq)];
        gradf = gradf + phi_bar.'*Bq1*phi_bar*cq1 + phi_bar.'*cq1*Bq1*phi_bar ...
            + cq1.'*phi_bar*Bq1.'*phi_bar + (phi_bar.'*Bq2*phi_bar*cq2 + phi_bar.'*cq2*Bq2*phi_bar ...
            + cq2.'*phi_bar*Bq2.'*phi_bar);
        Hessf = Hessf + cq1*phi_bar.'*Bq1.' + cq1*phi_bar.'*Bq1 + Bq1*phi_bar*cq1.' ...
            + phi_bar.'*cq1*Bq1 + Bq1.'*phi_bar*cq1.' + phi_bar.'*cq1*Bq1.' ...
            + (cq2*phi_bar.'*Bq2.' + cq2*phi_bar.'*Bq2 + Bq2*phi_bar*cq2.' ...
            + phi_bar.'*cq2*Bq2 + Bq2.'*phi_bar*cq2.' + phi_bar.'*cq2*Bq2.');
    end
    [~,eD] = eigsort(Hessf);
    lambda3 = max(real(diag(eD)));
    lt_bar = gradf-lambda3*phi_bar;
    ft = ft + U*ftv_bar + U*lt_bar-rho*z + mu2;
    Ft = Ft + (lambda2 + lambda3 + rho)/2*eye(N);
    for i = 1:1:N
        Ft(i,i) = real(Ft(i,i)) + 1e-6*real(Ft(1,1));
    end
    sc = norm(ft,2)^2;
    Ft = Ft/sc;
    ft = ft/sc;
    R = chol(Ft);
    %%%% solve phi
    cvx_begin quiet
    variable phi(N,1) complex
    minimize real(phi'*(R'*R)*phi) + real(phi'*ft)
    subject to
    real(gamma_CI(:,1).*(Hkl*x + Gkl_tilde*phi)) >= 1;
    real(gamma_CI(:,2).*(Hkl*x + Gkl_tilde*phi)) >= 1;
    abs(phi) <= 1;
    cvx_end

    Phi2 = phi*phi.';
    phi_bar = [real(phi);imag(phi)];
    Fqphi = zeros(M*L,M*L,Q+1);
    Fqphi(:,:,end) = A(:,:,end) + kron(eye(L),ht'*phi.'*B(:,:,end)) + kron(eye(L),B(:,:,end)'*phi*ht)...
        + kron(eye(L),B(:,:,end)'*Phi2*B(:,:,end));
    st = Fqphi(:,:,end)*x;
    Mt = sigmar2*eye(M*L);
    for q = 1:1:Q
        Fqphi(:,:,q) = ( A(:,:,q) + kron(eye(L),Hc(q,:)'*phi.'*B(:,:,q)) + kron(eye(L),B(:,:,q)'*phi*Hc(q,:))...
            + kron(eye(L),B(:,:,q)'*Phi2*B(:,:,q)) )*Jr(:,:,q);
        Mt = Mt + sigma2*Fqphi(:,:,q)*X*Fqphi(:,:,q)';
    end
    mst = Mt\st;
    Mst = mst*mst';

    %%%% update y and z
    y = sqrt(P/(M*L))*exp(1i*angle(rho*x+mu1));
    z = exp(1i*angle(rho*phi+mu2));
    %%% update dual variables
    mu1 = mu1 + rho*(x-y);
    mu2 = mu2 + rho*(phi-z);

    Dt = zeros(M*L,M*L);
    for q = 1:1:Q
        Dt = Dt + 2*sigma2*Fqphi(:,:,q)'*Mst*Fqphi(:,:,q);
    end
    Dt = Dt + rho/2*eye(M*L);
    dt = -2*Fqphi(:,:,end)'*mst-rho*y + mu1;

    Hkl_tilde = zeros(K*L,M*L);
    for ell = 1:1:L
        for k = 1:1:K
            Hkl_tilde((ell-1)*K+k,:) = Hkl((ell-1)*K+k,:) + kron(Iell(ell,:),phi.')*Gk(:,:,k);
        end
    end

    Vobj(iter) = real(-st'*mst) + 0.5*rho*(norm(phi-z+mu2/rho,2))^2 + 0.5*rho*(norm(x-y+mu1/rho,2))^2;
    VSINR(iter) = 10*log10(real(sigma2*st'*mst));
    Vres(iter) = norm(x-y,2)^2 + norm(phi-z,2)^2;
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


