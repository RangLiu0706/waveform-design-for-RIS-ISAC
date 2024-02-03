% This Matlab script can be used to generate Fig. 8 in the paper:
% R. Liu, M. Li, Y. Liu, Q. Wu, and Q. Liu, “Joint transmit waveform and passive beamforming design for RIS-aided DFRC systems,”IEEE J. Sel. Topics Signal Process., vol. 16, no .5, pp. 995-1010, Aug. 2022.
% Download this paper at: https://ieeexplore.ieee.org/document/9769997
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-01-28
clear
clc
%%% system settings
Prms.M = 6;  M = Prms.M; %%% number of transmit/receive antennas
Prms.N = 64; N = Prms.N; %%% number of RIS elements
Prms.K = 3;  K = Prms.K; %%% number of users
Prms.L = 20; L = Prms.L; %%% number of samples
Prms.clutter = [0 -50; 1 -10; 2 40]; clutter = Prms.clutter; %%% range-angle of clutter
Prms.Q = 3; Q = Prms.Q; %%% number of clutter patches
Prms.sigmar2 = 10^(-11); sigmar2 = Prms.sigmar2; %%% radar noise
Prms.sigmac2 = 10^(-11); sigmac2 = Prms.sigmac2; %%% communication noise
Prms.sigma2 = 1; sigma2 = Prms.sigma2; %%% RCS
Prms.Phi = pi/4; Phi = Prms.Phi;  %%% QPSK modulation
Prms.Nmax = 2000; Nmax = Prms.Nmax; %%% maximum iterations
Prms.res_th = 5e-4;%%% convergence tolerance
Prms.P = L*100; P = Prms.P; %%% total transmit power
%%% distances
dt = 30;
drt = 3;
dg = 30;
drk = 3;
dk = 30;
dcell = 50;
%%% path-loss
alpha_rt = 2.8;
alpha_rk = 2.8;
alpha_g = 2.5;

ht = exp(-1j*(0:1:M-1)*pi);
Channel.hrt = sqrt(10^(-3)*drt^(-alpha_rt))*exp(-1j*(0:1:N-1)*pi);
Hc = zeros(Q,M);
Channel.Hrc = zeros(Q,N);
for q = 1:1:Q
    theta = clutter(q,2);
    Hc(q,:) = exp(-1j*(0:1:M-1)*pi*sin(theta/180*pi));
    Channel.Hrc(q,:) = sqrt(10^(-3)*(drt+dcell*clutter(q,1))^(-alpha_rt))*exp(-1j*(0:1:N-1)*pi*sin(theta/180*pi)); %
end

SNR = 10*ones(1,K);
Prms.gamma = sqrt(sigmac2*10.^(0.1*SNR));

alpha_range = (2:0.4:4);

N_sim = 10000;
SINR_my_CI = zeros(1,length(alpha_range));
SINR_my_ZF = zeros(1,length(alpha_range));
SINR_my_MMSE = zeros(1,length(alpha_range));
SINR_radar_only = zeros(1,length(alpha_range));
SINR_wo_RIS_radar = zeros(1,length(alpha_range));
SINR_wo_RIS_CI = zeros(1,length(alpha_range));
SINR_wo_RIS_ZF = zeros(1,length(alpha_range));
SINR_wo_RIS_MMSE = zeros(1,length(alpha_range));

epsilon = 1e-9;

for sim = 1:1:N_sim
    tic
    sim
    Channel.G = sqrt(10^(-3)*dg^(-alpha_g))*exp(-1j*(0:1:N-1)'*pi*sin(pi*rand-pi/2))*exp(-1j*(0:1:M-1)*pi*sin(pi*rand-pi/2));
    Hu = (rand(K,M) + 1i*rand(K,M))/sqrt(2);%%sqrt(10^(-3)*dk^(-alpha_k))
    Channel.Hru = sqrt(10^(-3)*drk^(-alpha_rk))*(rand(K,N) + 1i*rand(K,N))/sqrt(2);
    S = exp(1i*(Phi+2*Phi*randi([0,pi/Phi-1],K,L)));

    phi = exp(1i*2*pi*rand(N,1));
    for alpha_index = 1:1:length(alpha_range)

        alpha_t = alpha_range(alpha_index)
        alpha_k = alpha_t;

        Channel.ht = ht*sqrt(10^(-3)*dt^(-alpha_t));
        for q = 1:1:Q
            Channel.Hc(q,:) = Hc(q,:)*sqrt(10^(-3)*(dt+dcell*clutter(q,1))^(-alpha_t));
        end
        Channel.Hu = Hu*sqrt(10^(-3)*dk^(-alpha_k));

        [x_radar,phi_radar,VSINR_radar] = get_x_RIS_radar(Prms,Channel);
        SINR_radar_only(alpha_index) = SINR_radar_only(alpha_index) + VSINR_radar(end);

        [x_radar_wo_RIS,VSINR_radar_wo_RIS] = get_x_woRIS_radar(Prms,Channel);
        SINR_wo_RIS_radar(alpha_index) = SINR_wo_RIS_radar(alpha_index) + VSINR_radar_wo_RIS(end);

        [x_my,phi_my,VSINR_my] = get_x_phi_CI(Prms,Channel,S);
        SINR_my_CI(alpha_index) = SINR_my_CI(alpha_index) + VSINR_my(end);

        [x_wo_RIS,VSINR_wo_RIS] = get_x_woRIS_CI(Prms,Channel,S);
        SINR_wo_RIS_CI(alpha_index) = SINR_wo_RIS_CI(alpha_index) + VSINR_wo_RIS(end);

        [x_ZF,phi_ZF,VSINR_ZF] = get_x_phi_ZF(Prms,Channel,S);
        SINR_my_ZF(alpha_index) = SINR_my_ZF(alpha_index) + VSINR_ZF(end);

        [x_wo_RIS_ZF,VSINR_wo_RIS_ZF] = get_x_woRIS_ZF(Prms,Channel,S);
        SINR_wo_RIS_ZF(alpha_index) = SINR_wo_RIS_ZF(alpha_index) + VSINR_wo_RIS_ZF(end);

        [x_MMSE,phi_MMSE,VSINR_MMSE] = get_x_phi_MMSE(Prms,Channel,S,epsilon);
        SINR_my_MMSE(alpha_index) = SINR_my_MMSE(alpha_index) + VSINR_MMSE(end);

        [x_wo_RIS_MMSE,VSINR_wo_RIS_MMSE] = get_x_woRIS_MMSE(Prms,Channel,S,epsilon);
        SINR_wo_RIS_MMSE(alpha_index) = SINR_wo_RIS_MMSE(alpha_index) + VSINR_wo_RIS_MMSE(end);

    end
    toc
end

SINR_my_CI = SINR_my_CI/sim;
SINR_my_ZF = SINR_my_ZF/sim;
SINR_my_MMSE = SINR_my_MMSE/sim;
SINR_radar_only = SINR_radar_only/sim;
SINR_wo_RIS_radar = SINR_wo_RIS_radar/sim;
SINR_wo_RIS_CI = SINR_wo_RIS_CI/sim;
SINR_wo_RIS_ZF = SINR_wo_RIS_ZF/sim;
SINR_wo_RIS_MMSE = SINR_wo_RIS_MMSE/sim;

figure
plot(alpha_range,SINR_my_CI,'-o','color',[0.05,0.4,0.05],'LineWidth',1.5)
hold on
plot(alpha_range,SINR_wo_RIS_CI,'--o','color',[0.05,0.4,0.05],'LineWidth',1.5)
plot(alpha_range,SINR_my_MMSE,'-s','color',[0.8,0,0],'LineWidth',1.5)
plot(alpha_range,SINR_wo_RIS_MMSE,'--s','color',[0.8,0,0],'LineWidth',1.5)
plot(alpha_range,SINR_my_ZF,'-^','color',[0,0,0.8],'LineWidth',1.5)
plot(alpha_range,SINR_wo_RIS_ZF,'--^','color',[0,0,0.8],'LineWidth',1.5)
plot(alpha_range,SINR_radar_only,'-x','color',[0,0,0],'LineWidth',1.5)
plot(alpha_range,SINR_wo_RIS_radar,'--x','color',[0,0,0],'LineWidth',1.5)
hold off
xlabel('Path-loss exponent {\it \alpha}_t');
ylabel('Radar SINR (dB)');
grid on
legend('CI-DFRC, w/ RIS','CI-DFRC, w/o RIS','MMSE-DFRC, w/ RIS','MMSE-DFRC, w/o RIS',...
    'ZF-DFRC, w/ RIS','ZF-DFRC, w/o RIS','Radar-only, w/ RIS', 'Radar-only, w/o RIS');
