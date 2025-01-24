function [ewma_stat, upper_limit] = EWMA_one_sided_limit(stats, lambda, L_para, miu_0, sigma_stat)
    % 计算任意统计量的单边EWMA控制上限
    %
    % 输入参数：
    % stats - 输入的统计量序列 (向量)
    % lambda - 平滑系数 (通常取值在0.05到0.25之间)
    % L_para - 控制限倍数
    % mu_0 - 统计量的均值
    % sigma_stat - 统计量的标准差
    %
    % 输出参数：
    % ewma_stat - 计算得到的EWMA统计量
    % upper_limit - 对应的上控制限

    testing_epoch = length(stats);  % 测试周期数
    upper_limit = zeros(testing_epoch, 1);
    ewma_stat = zeros(testing_epoch, 1);

    for i = 1:testing_epoch
        tmp = L_para * sigma_stat * sqrt(lambda / (2 - lambda) * (1 - (1 - lambda)^(2 * i)));
        upper_limit(i) = miu_0 + 2 * tmp;

        if i == 1
            ewma_stat(i) = lambda * stats(i) + (1 - lambda) * miu_0;
        else
            ewma_stat(i) = lambda * stats(i) + (1 - lambda) * ewma_stat(i-1);
        end
    end
end
