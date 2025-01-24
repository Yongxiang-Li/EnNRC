function [ X, Qyt ] = mul_ennrc( signal, L, Periods )
%MUL_ENNRC calculate the statistics and control limit
%  Consider the correlation
    signal = zscore(signal);
    Qyt = zeros(length(Periods), 1);
    for i = 1:length(Periods)
        p = Periods(i);
        Qyt(i) = en_nrc(signal, p);
    end

    Z = L \ Qyt;
    Qyt = flipud(Qyt);
    X = sum(Z.^2)'; % chi square distribution with degree of freedom length(period)

end

