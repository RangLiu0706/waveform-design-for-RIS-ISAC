% This Matlab script can be used to generate Fig. 6 in the paper:
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
%%%% channel settings
%%% distances
dt = 30;
drt = 3;
dg = 30;
drk = 3;
dk = 30;
dcell = 50;
%%% path-loss
alpha_t = 3;
alpha_rt = 2.8;
alpha_k = 3;
alpha_rk = 2.8;
alpha_g = 2.5;

Channel.ht = sqrt(10^(-3)*dt^(-alpha_t))*exp(-1j*(0:1:M-1)*pi);
Channel.hrt = sqrt(10^(-3)*drt^(-alpha_rt))*exp(-1j*(0:1:N-1)*pi);
Channel.Hc = zeros(Q,M);
Channel.Hrc = zeros(Q,N);
for q = 1:1:Q
    theta = clutter(q,2);
    Channel.Hc(q,:) = sqrt(10^(-3)*(dt+dcell*clutter(q,1))^(-alpha_t))*exp(-1j*(0:1:M-1)*pi*sin(theta/180*pi));%/sqrt(M); %
    Channel.Hrc(q,:) = sqrt(10^(-3)*(drt+dcell*clutter(q,1))^(-alpha_rt))*exp(-1j*(0:1:N-1)*pi*sin(theta/180*pi)); %
end

SNR_dB_range = (10:5:35);
N_sim = 10000;
SINR_my_CI = zeros(1,length(SNR_dB_range));
SINR_my_ZF = zeros(1,length(SNR_dB_range));
SINR_my_MMSE = zeros(1,length(SNR_dB_range));
SINR_RIS_radar = zeros(1,length(SNR_dB_range));
SINR_wo_RIS_radar = zeros(1,length(SNR_dB_range));
SINR_wo_RIS_DFRC = zeros(1,length(SNR_dB_range));
SINR_wo_RIS_ZF_DFRC = zeros(1,length(SNR_dB_range));
SINR_wo_RIS_MMSE_DFRC = zeros(1,length(SNR_dB_range));

epsilon = 1e-9;

for sim = 1:1:N_sim
    tic
    sim
    Channel.G = sqrt(10^(-3)*dg^(-alpha_g))*exp(-1j*(0:1:N-1)'*pi*sin(pi*rand-pi/2))*exp(-1j*(0:1:M-1)*pi*sin(pi*rand-pi/2));%/sqrt(M);
    Channel.Hu = sqrt(10^(-3)*dk^(-alpha_k))*(rand(K,M) + 1i*rand(K,M))/sqrt(2);%/sqrt(M);
    Channel.Hru = sqrt(10^(-3)*drk^(-alpha_rk))*(rand(K,N) + 1i*rand(K,N))/sqrt(2);
    S = exp(1i*(Phi+2*Phi*randi([0,pi/Phi-1],K,L)));

    [x_radar,phi_radar,VSINR_radar] = get_x_RIS_radar(Prms,Channel);
    SINR_RIS_radar = SINR_RIS_radar + VSINR_radar(end)*ones(1,length(SNR_dB_range));

    [x_radar_wo_RIS,VSINR_radar_wo_RIS] = get_x_woRIS_radar(Prms,Channel);
    SINR_wo_RIS_radar = SINR_wo_RIS_radar + VSINR_radar_wo_RIS(end)*ones(1,length(SNR_dB_range));

    phi = exp(1i*2*pi*rand(N,1));
    for SNR_index = 1:1:length(SNR_dB_range)
        SNR = SNR_dB_range(SNR_index)*ones(1,K)
        Prms.gamma = sqrt(sigmac2*10.^(0.1*SNR));

        [x_my,phi_my,VSINR_my] = get_x_phi_CI(Prms,Channel,S);
        SINR_my_CI(SNR_index) = SINR_my_CI(SNR_index) + VSINR_my(end);

        [x_wo_RIS,VSINR_wo_RIS] = get_x_woRIS_CI(Prms,Channel,S);
        SINR_wo_RIS_DFRC(SNR_index) = SINR_wo_RIS_DFRC(SNR_index) + VSINR_wo_RIS(end);

        [x_ZF,phi_ZF,VSINR_ZF] = get_x_phi_ZF(Prms,Channel,S);
        SINR_my_ZF(SNR_index) = SINR_my_ZF(SNR_index) + VSINR_ZF(end);

        [x_wo_RIS_ZF,VSINR_wo_RIS_ZF] = get_x_woRIS_ZF(Prms,Channel,S);
        SINR_wo_RIS_ZF_DFRC(SNR_index) = SINR_wo_RIS_ZF_DFRC(SNR_index) + VSINR_wo_RIS_ZF(end);

        [x_MMSE,phi_MMSE,VSINR_MMSE] = get_x_phi_MMSE(Prms,Channel,S,epsilon);
        SINR_my_MMSE(SNR_index) = SINR_my_MMSE(SNR_index) + VSINR_MMSE(end);

        [x_wo_RIS_MMSE,VSINR_wo_RIS_MMSE] = get_x_woRIS_MMSE(Prms,Channel,S,epsilon);
        SINR_wo_RIS_MMSE_DFRC(SNR_index) = SINR_wo_RIS_MMSE_DFRC(SNR_index) + VSINR_wo_RIS_MMSE(end);

    end
    toc
end

SINR_my_CI = SINR_my_CI/sim;
SINR_my_ZF = SINR_my_ZF/sim;
SINR_my_MMSE = SINR_my_MMSE/sim;
SINR_RIS_radar = SINR_RIS_radar/sim;
SINR_wo_RIS_radar = SINR_wo_RIS_radar/sim;
SINR_wo_RIS_DFRC = SINR_wo_RIS_DFRC/sim;
SINR_wo_RIS_ZF_DFRC = SINR_wo_RIS_ZF_DFRC/sim;
SINR_wo_RIS_MMSE_DFRC = SINR_wo_RIS_MMSE_DFRC/sim;

figure
plot(SNR_dB_range,SINR_my_CI,'-o','color',[0.05,0.4,0.05],'LineWidth',1.5)
hold on
plot(SNR_dB_range,SINR_wo_RIS_DFRC,'--o','color',[0.05,0.4,0.05],'LineWidth',1.5)
plot(SNR_dB_range,SINR_my_MMSE,'-s','color',[0.8,0,0],'LineWidth',1.5)
plot(SNR_dB_range,SINR_wo_RIS_MMSE_DFRC,'--s','color',[0.8,0,0],'LineWidth',1.5)
plot(SNR_dB_range,SINR_my_ZF,'-^','color',[0,0,0.8],'LineWidth',1.5)
plot(SNR_dB_range,SINR_wo_RIS_ZF_DFRC,'--^','color',[0,0,0.8],'LineWidth',1.5)
plot(SNR_dB_range,SINR_RIS_radar,'-x','color',[0,0,0],'LineWidth',1.5)
plot(SNR_dB_range,SINR_wo_RIS_radar,'--x','color',[0,0,0],'LineWidth',1.5)
hold off
xlabel('Communication QoS \Gamma (dB)');
ylabel('Radar SINR (dB)');
grid on
legend('CI-DFRC, w/ RIS','CI-DFRC, w/o RIS','MMSE-DFRC, w/ RIS','MMSE-DFRC, w/o RIS',...
    'ZF-DFRC, w/ RIS','ZF-DFRC, w/o RIS','Radar-only, w/ RIS', 'Radar-only, w/o RIS');

