function [miu_0, sigma_stat] = cal_miu_sig(num_of_TE, len_of_signal, train_signal, L, Periods)
%CAL_MIU_SIG Estimate mu and sigma
%   此处显示详细说明


stat = zeros(num_of_TE, 1);
for i = 0:(num_of_TE - 1)
    startindex = 1 + len_of_signal*i;
    endindex = len_of_signal*i + len_of_signal;
    signal = train_signal(startindex:endindex);
    [stat(i+1),~] = Max_NRC(signal, L, Periods);
end

miu_0 = mean(stat);
sigma_stat = std(stat);
end

