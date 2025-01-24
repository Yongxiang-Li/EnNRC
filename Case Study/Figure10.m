close all;
folderPath = '../Data/Bearing1_3';
addpath('../routine'); 
addpath('../Data/Bearing1_3'); % Data are publicly available at https://biaowang.tech/xjtu-sy-bearing-datasets
rng('default')

sigma = 1; 
lambda = 0.25; 
L_para = 3; % L_para determines the limit rule: L_para=3 means 3-sigma
Periods = 3277:-1:2;
P = length(Periods);
len = 32768;

num = 158;  % number of epochs
mod = 0;
training_epoch = 20;
testing_epoch = num - mod;

% Searching range for MO statistic
if ~exist('C', 'var')
    C = sigma^4*get_cov_matrix(len, Periods);
    L = chol(C, 'lower');
end

% Initiate Statistics
t2stat = zeros(training_epoch, 1);
maxstat = zeros(training_epoch, 1);
fstat = zeros(training_epoch, 1);

for i = 1:training_epoch
    filename = [num2str(i + mod) '.csv'];
    
    [data] = csvread(filename,1,0);
    signal = data(:,1).^2 + data(:,2).^2;
    signal = zscore(signal);
    
    [t2stat(i), ~] = mul_ennrc(signal, L, Periods);
    [maxstat(i),~,~] = Max_NRC(signal, L, Periods);
    [fstat(i),~] = Fchart(signal);
end

% sample mean and std
mu_t2 = mean(t2stat);
sigma_t2 = std(t2stat);
mu_max = mean(maxstat);

sigma_max = std(maxstat);
mu_f = mean(fstat);
sigma_f = std(fstat);

test_stat_t2 = 0;
test_stat_max = 0;
test_stat_f = 0;

% Phase 2: Calculate the limit and EWMA statistic simultaneously
upper_limit_t2 = zeros(testing_epoch, 1);
ewma_stat_t2 = zeros(testing_epoch, 1);
upper_limit_max = zeros(testing_epoch, 1);
ewma_stat_max = zeros(testing_epoch, 1);
upper_limit_f = zeros(testing_epoch, 1);
ewma_stat_f = zeros(testing_epoch, 1);

for i = 1:testing_epoch
    filename = [num2str(i + mod) '.csv'];
    [data] = csvread(filename,1,0);
    signal = data(:,1).^2 + data(:,2).^2;
    signal = zscore(signal);
    
    % control limit
    tmp_t2 = L_para * sigma_t2 * sqrt(lambda/(2-lambda)*(1-(1-lambda)^(2*i)));
    upper_limit_t2(i) = mu_t2 + 2 * tmp_t2;
    
    tmp_max = L_para * sigma_max * sqrt(lambda/(2-lambda)*(1-(1-lambda)^(2*i)));
    upper_limit_max(i) = mu_max + 2 * tmp_max;
    
    tmp_f = L_para * sigma_f * sqrt(lambda/(2-lambda)*(1-(1-lambda)^(2*i)));
    upper_limit_f(i) = mu_f + 2 * tmp_f;
    
    [test_stat_t2, ~] = mul_ennrc(signal, L, Periods);
    [test_stat_max, ~, ~] = Max_NRC(signal, L, Periods);
    [test_stat_f, ~] = Fchart(signal);
    
    % compute EWMA statistics
    if i == 1
        ewma_stat_t2(i) = lambda * test_stat_t2 + (1 - lambda) * mu_t2;
        ewma_stat_max(i) = lambda * test_stat_max + (1 - lambda) * mu_max;
        ewma_stat_f(i) = lambda * test_stat_f + (1 - lambda) * mu_f;
    else
        ewma_stat_t2(i) = lambda * test_stat_t2 + (1 - lambda) * ewma_stat_t2(i-1);
        ewma_stat_max(i) = lambda * test_stat_max + (1 - lambda) * ewma_stat_max(i-1);
        ewma_stat_f(i) = lambda * test_stat_f + (1 - lambda) * ewma_stat_f(i-1);
    end
end

% plot
figure;
plot((1:testing_epoch)', ewma_stat_max,'Color',[0.85,0.33,0.10]);
hold on;
plot((1:testing_epoch)', upper_limit_max, '-','color', [0.00  0.45  0.74]);
hold on
xlabel('Running Epoch');
ylabel('MO EWMA Control Chart');


