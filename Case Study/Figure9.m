close all;
addpath('../routine/'); 
addpath('../Data/Bearing1_3'); % Data are publicly available at https://biaowang.tech/xjtu-sy-bearing-datasets
rng('default')

sigma = 1;    Periods = 3277:-1:2;
P = length(Periods);
num = 158;
t2stat = zeros(num,1);    maxstat = zeros(num,1);    fstat = zeros(num,1);

for i = 1:num
    filename = [num2str(i) '.csv'];
    [data]= csvread(filename,1,0);

    signal = data(:,1).^2 + data(:,2).^2;
    signal = zscore(signal);
    
    len = length(signal);
    if ~exist('C', 'var')
        C = sigma^4*get_cov_matrix(len, Periods);
        L = chol(C, 'lower');
    end
    
    [t2stat(i), ~] = mul_ennrc(signal, L, Periods);
    [maxstat(i),~,~] = Max_NRC(signal, L, Periods);
    [fstat(i),~] = Fchart(signal);
end

% Calculate the limit
pr_limit = 0.0001;    len = length(signal);
upper_nrc = chi2inv(1 - pr_limit, P);
upper_max = MaxNRC_CL(P, pr_limit);
upper_fisher = Fchart_CL(floor(len/2),pr_limit);

% Plot the result
figure; 
subplot(1,3,1)
plot(t2stat)
line([0,num],[upper_nrc,upper_nrc],'linestyle','--','Color',[1,0,0]);
xlabel('Running Epoch')
ylabel('$\textrm{T}^2$ Control Chart','Interpreter','latex')
xlim([0 num])
RemoveSubplotWhiteArea(gca, 1, 3, 1, 1, 0.05);

subplot(1,3,2)
plot(maxstat)
line([0,num],[upper_max,upper_max],'linestyle','--','Color',[1,0,0]);
xlabel('Running Epoch')
ylabel('MO Control Chart')
xlim([0 num])
RemoveSubplotWhiteArea(gca, 1, 3, 1, 2, 0.05);

subplot(1,3,3)
plot(fstat)
line([0,num],[upper_fisher,upper_fisher],'linestyle','--','Color',[1,0,0]);
xlabel('Running Epoch')
ylabel('Spectral Control Chart')
RemoveSubplotWhiteArea(gca, 1, 3, 1, 3, 0.05);