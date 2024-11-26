close all;
addpath('../routine')
rng('default')
fs = 100;    f0 = 40;    T = 1;
zeta0 = 0.01;    sigma = 1; % if no noise                
pr_limit = 0.001;    
len = 5000;      SNR = 0;     step = 100;

P = floor(len/10);
Periods = (P+1):-1:2;
C = get_cov_matrix(len, Periods);
L = chol(C, 'lower');

White_signal = randn(len*200,1);

GP_signal = sigma*randGP(len*200, 0.3); % get signal
noise = randn(size(GP_signal));
noise = noise*rms(GP_signal)/rms(noise)/10^(SNR/20); % 20*log10(rms(originSignal)/rms(noise))
GP_signal = GP_signal + noise;

SNR = -13;
[Periodic_signal, t] = get_normal_transient_signal(len*200, f0, zeta0, T, fs, 0); % get signal
noise = randn(size(Periodic_signal));
noise = noise*rms(Periodic_signal)/rms(noise)/10^(SNR/20); % 20*log10(rms(originSignal)/rms(noise))
Periodic_signal = Periodic_signal + noise;

originSignal = [White_signal; GP_signal; Periodic_signal];

upper_nrc = chi2inv(1 - pr_limit, P);
% upper_mcnrc = NRC_mcCL(len, Periods, pr_limit, L);

upper_max = MaxNRC_CL(P, pr_limit);
% upper_mcmax = MaxNRC_MCCL( len, L, Periods, pr_limit );

upper_fisher = Fchart_CL(floor(len/2),pr_limit);

tic
for i = 0:floor((length(originSignal)-len)/len)
    startindex = 1 + len*i;
    endindex = len*i + len;
    signal = originSignal(startindex:endindex);
    [nrcstat(i+1), ~] = mul_ennrc(signal, L, Periods);
    [maxstat(i+1),~] = Max_NRC(signal, L, Periods);
    fstat(i+1) = Fchart(signal);
end
toc

nrcstat = nrcstat(:);
maxstat = maxstat(:);
fstat = fstat(:);

figure;
subplot(1,3,1)
plot((1:200)',nrcstat(1:200),'Color',[0.85,0.33,0.10]);
hold on
plot((200:400)',nrcstat(200:400),'Color',[0.93,0.69,0.13])
hold on
plot((400:600)',nrcstat(400:600),'Color',[0.72,0.27,1])
hold on
plot(upper_nrc*ones(size(nrcstat)), '-','color', [0.00  0.45  0.74])
hold on
xlabel('Running Epoch')
ylabel('$\textrm{T}^2$ Control Chart','Interpreter','latex')
xlim([0 600])
RemoveSubplotWhiteArea(gca, 1, 3, 1, 1, 0.05);

subplot(1,3,2)
plot((1:200)',maxstat(1:200),'Color',[0.85,0.33,0.10]);
hold on
plot((200:400)',maxstat(200:400),'Color',[0.93,0.69,0.13])
hold on
plot((400:600)',maxstat(400:600),'Color',[0.72,0.27,1])
hold on
plot(upper_max*ones(size(maxstat)), '-','color', [0.00  0.45  0.74])
hold on
xlabel('Running Epoch')
ylabel('MO Control Chart')
xlim([0 600])
RemoveSubplotWhiteArea(gca, 1, 3, 1, 2, 0.05);

subplot(1,3,3)
plot((1:200)',fstat(1:200),'Color',[0.85,0.33,0.10]);
hold on
plot((200:400)',fstat(200:400),'Color',[0.93,0.69,0.13])
hold on
plot((400:600)',fstat(400:600),'Color',[0.72,0.27,1])
hold on
plot(upper_fisher*ones(size(fstat)), '-','color', [0.00  0.45  0.74])
hold on
xlabel('Running Epoch')
ylabel('Fisher Control Chart')
xlim([0 600])
RemoveSubplotWhiteArea(gca, 1, 3, 1, 3, 0.05);


