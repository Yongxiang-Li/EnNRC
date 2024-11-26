function [ control_limit ] = MaxNRC_MCCL( len, L, Periods, pr )
%   Calculate the control limit for the Max NRC statistic
%   Suppose pr is small, which is ususally the case
    pr = 1 - pr;
    X = zeros(10000,1);
    for i = 1:10000
        signal = randn(len,1);
        [X(i),~,~] = Max_NRC( signal, L, Periods );
    end
    control_limit = quantile(X,pr);
end

