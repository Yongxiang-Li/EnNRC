function [miu_0, sigma_stat] = cal_miu_sig_from_stats(stats)
%CAL_MIU_SIG_FROM_STATS Calculate the mean (mu_0) and standard deviation (sigma_stat) from a set of statistics
%   该函数接收一组统计量并计算其均值和标准差

% 输入参数：
% stats - 输入的统计量序列 (向量)

% 输出参数：
% miu_0 - 统计量的均值
% sigma_stat - 统计量的标准差

% 计算均值和标准差
train_epoch = length(stats) / 10;
train_stats = stats(1: train_epoch);
miu_0 = mean(train_stats);
sigma_stat = std(train_stats);
end